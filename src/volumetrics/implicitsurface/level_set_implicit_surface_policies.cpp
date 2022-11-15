// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/math/reconstruction_filters.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>
#include <frantic/volumetrics/implicitsurface/level_set_implicit_surface_policies.hpp>
#include <frantic/volumetrics/levelset/rle_trilerp.hpp>

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using namespace frantic::geometry;
using namespace frantic::volumetrics::levelset;

namespace frantic {
namespace volumetrics {
namespace implicitsurface {

direct_linear_rle_level_set_is_policy::direct_linear_rle_level_set_is_policy( const levelset::rle_level_set& ls,
                                                                              int boundsClip )
    : m_levelSet( ls )
    , m_boundsClip( boundsClip ) {}

const voxel_coord_system& direct_linear_rle_level_set_is_policy::get_voxel_coord_system() const {
    return m_levelSet.get_voxel_coord_system();
}

// This returns the XY bounds within which to operate
boundbox3 direct_linear_rle_level_set_is_policy::get_voxel_bounds() const {
    boundbox3 result = m_levelSet.get_rle_index_spec().outer_bounds();
    if( !m_boundsClip )
        result.expand( 3 );
    else
        result.expand( -m_boundsClip + 3 );
    return result;
}

// For the purposes of marching cubes, the samples of the level set are on the voxel corners, so the voxels are
// shifted by 1/2 from the alternate centered interpretation.

void direct_linear_rle_level_set_is_policy::fill_voxel_corner_densities(
    const boundrect2& xyExtents, int z, std::vector<float>& outVoxelCornerValues ) const {
    m_levelSet.fill_plane( xyExtents, z, outVoxelCornerValues );
}

void direct_linear_rle_level_set_is_policy::fill_sparse_voxel_corner_densities(
    const frantic::graphics2d::boundrect2& xyExtents, int z, float* outVoxelCornerValues,
    frantic::volumetrics::rle_plane& rlp ) {

    std::vector<frantic::tstring> channelNames;
    channelNames.push_back( _T("SignedDistance") );
    std::vector<char*> channelData;
    channelData.push_back( (char*)outVoxelCornerValues );
    m_levelSet.fill_sparse_plane_channel_data( xyExtents, z, channelNames, channelData, rlp );
    // m_levelSet.fill_sparse_plane( xyExtents, z, outVoxelCornerValues, rlp );
}

void direct_linear_rle_level_set_is_policy::fill_sparse_voxel_corner_densities(
    const frantic::graphics2d::boundrect2& xyExtents, int z, float* outVoxelCornerValues, detail::empty_workspace& ) {
    m_levelSet.fill_plane( xyExtents, z, outVoxelCornerValues );
}

void direct_linear_rle_level_set_is_policy::fill_sparse_channel_data( const frantic::graphics2d::boundrect2& xyExtents,
                                                                      int z,
                                                                      std::vector<frantic::tstring>& channelNames,
                                                                      std::vector<char*>& channelData,
                                                                      frantic::volumetrics::rle_plane& rlp ) {

    m_levelSet.fill_sparse_plane_channel_data( xyExtents, z, channelNames, channelData, rlp );
}

vector3f direct_linear_rle_level_set_is_policy::corner_sample_coord_to_world( const vector3& voxelCornerCoord ) const {
    return m_levelSet.get_voxel_coord_system().get_world_voxel_center( voxelCornerCoord );
}

/**
 * This function finds the location of the isosurface between the two voxel corners provided, and sets the internal
 * state of this is policy based on that.
 *
 * @note More advanced versions will do a numerical solve between the world-space
 *       coordinate versions of voxelCorner0 and voxelCorner1
 *
 * @param  density0  The density function at voxelCorner0.
 * @param  density1  The density function at voxelCorner1.
 * @param  voxelCorner0  The first of two corners between which the isosurface is contained.
 * @param  voxelCorner1  The second of two corners between which the isosurface is contained.
 */
void direct_linear_rle_level_set_is_policy::find_edge_isosurface_location( float density0, float density1,
                                                                           const vector3& voxelCorner0,
                                                                           const vector3& voxelCorner1 ) {
    float densityDelta = density1 - density0;
    if( densityDelta == 0 )
        m_isosurfaceLocationAlpha = 0.5f;
    else
        m_isosurfaceLocationAlpha = fabsf( density0 / densityDelta );
    m_isosurfaceLocationVoxelCorner0 = voxelCorner0;
    m_isosurfaceLocationVoxelCorner1 = voxelCorner1;
    m_isosurfaceLocationWorld = ( 1 - m_isosurfaceLocationAlpha ) * corner_sample_coord_to_world( voxelCorner0 ) +
                                m_isosurfaceLocationAlpha * corner_sample_coord_to_world( voxelCorner1 );
}

/**
 * This function computes an additional vertex location internal to a voxel given two
 * opposing voxel corners.  It doesn't modify the internal state of the is policy.
 *
 * @param  density0  The density function at voxelCorner0.
 * @param  density1  The density function at voxelCorner1.
 * @param  voxelCorner0  The first of two corners between which the isosurface is contained.
 * @param  voxelCorner1  The second of two corners between which the isosurface is contained.
 * @return vector3f	 The 3D location of the vertex on the isosurface internal to the voxel.
 */
vector3f direct_linear_rle_level_set_is_policy::find_internal_isosurface_location( float density0, float density1,
                                                                                   const vector3& voxelCorner0,
                                                                                   const vector3& voxelCorner1 ) {

    float densityDelta = density1 - density0;
    float isosurfaceLocationAlpha;

    if( densityDelta == 0 )
        isosurfaceLocationAlpha = 0.5f;
    else
        isosurfaceLocationAlpha = fabsf( density0 / densityDelta );
    vector3 isosurfaceLocationVoxelCorner0 = voxelCorner0;
    vector3 isosurfaceLocationVoxelCorner1 = voxelCorner1;
    return ( 1 - isosurfaceLocationAlpha ) * corner_sample_coord_to_world( voxelCorner0 ) +
           isosurfaceLocationAlpha * corner_sample_coord_to_world( voxelCorner1 );
}

const vector3f& direct_linear_rle_level_set_is_policy::get_isosurface_location_world() const {
    return m_isosurfaceLocationWorld;
}

////////
// Functions for dealing with the named channels
////////

/**
 * This retrieves the named channels in the input data structure.
 *
 * @param  outNames  A vector of strings into which the channel names are added.
 */
void direct_linear_rle_level_set_is_policy::get_channel_names( std::vector<frantic::tstring>& outNames ) const {
    m_levelSet.get_channel_names( outNames );
}

/**
 *	This retrieves the named channels in the level set that the policy is operating on.
 *	It adds the names to the existing vector, so if you want just the names this function provides,
 *	make sure to initialize the output vector to empty yourself.  It uses the channel propagation
 *	policy to determine which channels should be returned.
 *
 *	@param	outNames	A vector of channel name strings.
 */
void direct_linear_rle_level_set_is_policy::get_channel_names( const frantic::channels::channel_propagation_policy& cpp,
                                                               std::vector<frantic::tstring>& outNames ) const {
    std::vector<frantic::tstring> channelNames;
    m_levelSet.get_channel_names( channelNames );
    size_t channelCount = channelNames.size();

    // if ( cpp.is_channel_included("SignedDistance") ) {
    //	outNames.push_back("SignedDistance");
    // }

    for( size_t i = 0; i < channelCount; ++i ) {
        if( cpp.is_channel_included( channelNames[i] ) ) {
            outNames.push_back( channelNames[i] );
        }
    }
}

/**
 * This retrieves all the channel information associated with the named channels
 * in isp.
 *
 * @param	cpp						A channel propagation policy.
 * @param	outName					A vector of strings into which the channel names are added.
 * @param	outChannelType			A vector of data_type_t for each channel
 * @param	outChannelArity			A vector of size_t arities for each channel
 * @param	outChannelPrimitiveSize	A vector of size_t primitive sizes for each channel.
 */
void direct_linear_rle_level_set_is_policy::get_channel_info(
    frantic::channels::channel_propagation_policy cpp, std::vector<frantic::tstring>& outName,
    std::vector<frantic::channels::data_type_t>& outChannelType, std::vector<size_t>& outChannelArity,
    std::vector<size_t>& outChannelPrimitiveSize ) const {

    std::vector<frantic::tstring> channelNames;
    m_levelSet.get_channel_names( channelNames );
    size_t channelCount = channelNames.size();

    if( cpp.is_channel_included( _T("SignedDistance") ) ) {
        outName.push_back( _T("SignedDistance") );
        outChannelType.push_back( frantic::channels::data_type_float32 );
        outChannelArity.push_back( 1 );
        outChannelPrimitiveSize.push_back( 4 );
    }

    for( size_t i = 0; i < channelCount; ++i ) {
        if( cpp.is_channel_included( channelNames[i] ) ) {
            const_rle_channel_general_accessor acc = m_levelSet.get_channel_general_accessor( channelNames[i] );
            outName.push_back( channelNames[i] );
            outChannelType.push_back( acc.data_type() );
            outChannelArity.push_back( acc.arity() );
            outChannelPrimitiveSize.push_back( acc.primitive_size() );
        }
    }
}

/**
 * This function is called so that the channel accessors can be prepared,
 * and to report the data types for the algorithm to create corresponding
 * channels in the output data structure.
 *
 * @param  names  The channel names for which to prepare accessors.
 * @param  outDataTypes  The data types that the source data structure has for the corresponding names.
 */
void direct_linear_rle_level_set_is_policy::prepare_channel_accessors(
    const std::vector<frantic::tstring>& names, std::vector<std::pair<std::size_t, data_type_t>>& outDataTypes ) {
    outDataTypes.clear();
    outDataTypes.reserve( names.size() );
    m_inputChannels.clear();
    m_inputChannels.reserve( names.size() );
    for( unsigned i = 0; i < names.size(); ++i ) {
        m_inputChannels.push_back( m_levelSet.get_channel_general_accessor( names[i] ) );
        outDataTypes.push_back( std::make_pair( m_inputChannels.back().arity(), m_inputChannels.back().data_type() ) );
    }
}

/**
 * Given a set of writable trimesh3 accessors, this should add one data element to each of the channels,
 * which corresponds the the channel names given in prepare_channel_accessors.  The function uses the internal state of
 * the is policy for the location at which to look up the data
 *
 * @param  outputChannels  An array of trimesh3 vertex channel accessors.  One data element gets added to each.
 */
void direct_linear_rle_level_set_is_policy::add_vertex_to_channels(
    std::vector<geometry::trimesh3_vertex_channel_general_accessor>& outputChannels ) {
    // The same data indices are shared across all the channels, so we only have to get them once.
    // Special case when alpha is 0 or 1 to only look up one data entry
    int dataIndex0 = -1, dataIndex1 = -1;
    if( m_isosurfaceLocationAlpha > 0 )
        dataIndex0 = m_levelSet.get_rle_index_spec().XYZtoDataIndex( m_isosurfaceLocationVoxelCorner0 );
    if( m_isosurfaceLocationAlpha < 1 )
        dataIndex1 = m_levelSet.get_rle_index_spec().XYZtoDataIndex( m_isosurfaceLocationVoxelCorner1 );

    if( dataIndex0 >= 0 && dataIndex1 >= 0 ) {
        // Both data values are available, so do a linear interpolation
        for( unsigned i = 0; i < outputChannels.size(); ++i )
            m_inputChannels[i].get_linear_interpolated( dataIndex0, dataIndex1, m_isosurfaceLocationAlpha,
                                                        outputChannels[i].add_vertex() );
    } else if( dataIndex0 >= 0 ) {
        // Only the 0 data value is available, so just copy that value
        for( unsigned i = 0; i < outputChannels.size(); ++i )
            memcpy( outputChannels[i].add_vertex(), m_inputChannels[i].data( dataIndex0 ),
                    m_inputChannels[i].primitive_size() );
    } else if( dataIndex1 >= 0 ) {
        // Only the 1 data value is available, so just copy that value
        for( unsigned i = 0; i < outputChannels.size(); ++i )
            memcpy( outputChannels[i].add_vertex(), m_inputChannels[i].data( dataIndex1 ),
                    m_inputChannels[i].primitive_size() );
    } else {
        // There's no data available, so just add the vertex and set it to zero.
        for( unsigned i = 0; i < outputChannels.size(); ++i )
            memset( outputChannels[i].add_vertex(), 0, outputChannels[i].primitive_size() );
    }

    /*
    // DEBUGGING
    for(unsigned i = 0; i < outputChannels.size(); ++i ) {
      if( m_inputChannels[i].data_type() == channels::data_type_float32 &&
        !_finite(*reinterpret_cast<float*>(outputChannels[i].data(outputChannels[i].size()-1))) )
      {
        ofstream fout( "c:\\temp\\dbg_preview_linear.txt", ios::app|ios::out );
        fout << "Vertex " << outputChannels[i].size()-1 << " is bad\n";
        fout << "Alpha: " << m_isosurfaceLocationAlpha << "\n";
        fout << "Data indexes: " << dataIndex0 << " " << dataIndex1 << "\n";
        fout << "\n";
      }
    }
    */
}

void direct_linear_rle_level_set_is_policy::get_interface_widths( float& interfaceWidthInside,
                                                                  float& interfaceWidthOutside ) const {
    interfaceWidthInside = m_levelSet.get_interface_voxel_width_inside();
    interfaceWidthOutside = m_levelSet.get_interface_voxel_width_outside();
}

int direct_linear_rle_level_set_is_policy::get_exterior_region_code() const {
    return m_levelSet.get_rle_index_spec().get_exterior_region_code();
}

void direct_linear_rle_level_set_is_policy::add_edge_isosurface_vertex_to_mesh(
    float density0, float density1, const frantic::graphics::vector3& voxelCorner0,
    const frantic::graphics::vector3& voxelCorner1, std::size_t vertIndex,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels,
    frantic::geometry::trimesh3& outMesh, detail::empty_workspace& ) const {
    const float densityDelta = density1 - density0;

    float isosurfaceLocationAlpha;
    if( densityDelta == 0 )
        isosurfaceLocationAlpha = 0.5f;
    else
        isosurfaceLocationAlpha = fabsf( density0 / densityDelta );

    const vector3f isosurfaceLocationWorld =
        ( 1 - isosurfaceLocationAlpha ) * corner_sample_coord_to_world( voxelCorner0 ) +
        isosurfaceLocationAlpha * corner_sample_coord_to_world( voxelCorner1 );

    {
        // add a vert to the mesh at the surface location
        // tbb::spin_rw_mutex::scoped_lock meshLock(meshMutex, false);
        outMesh.vertices_ref()[vertIndex] = isosurfaceLocationWorld;
    }

    if( outputChannels.size() > 0 ) {
        // The same data indices are shared across all the channels, so we only have to get them once.
        // Special case when alpha is 0 or 1 to only look up one data entry
        int dataIndex0 = -1, dataIndex1 = -1;
        if( isosurfaceLocationAlpha > 0 )
            dataIndex0 = m_levelSet.get_rle_index_spec().XYZtoDataIndex( voxelCorner0 );
        if( isosurfaceLocationAlpha < 1 )
            dataIndex1 = m_levelSet.get_rle_index_spec().XYZtoDataIndex( voxelCorner1 );

        if( dataIndex0 >= 0 && dataIndex1 >= 0 ) {
            // tbb::spin_rw_mutex::scoped_lock meshLock(meshMutex, false);
            //  Both data values are available, so do a linear interpolation
            for( unsigned i = 0; i < outputChannels.size(); ++i )
                m_inputChannels[i].get_linear_interpolated( dataIndex0, dataIndex1, isosurfaceLocationAlpha,
                                                            outputChannels[i].data( vertIndex ) );
        } else if( dataIndex0 >= 0 ) {
            // tbb::spin_rw_mutex::scoped_lock meshLock(meshMutex, false);
            //  Only the 0 data value is available, so just copy that value
            for( unsigned i = 0; i < outputChannels.size(); ++i )
                memcpy( outputChannels[i].data( vertIndex ), m_inputChannels[i].data( dataIndex0 ),
                        m_inputChannels[i].primitive_size() );
        } else if( dataIndex1 >= 0 ) {
            // tbb::spin_rw_mutex::scoped_lock meshLock(meshMutex, false);
            //  Only the 1 data value is available, so just copy that value
            for( unsigned i = 0; i < outputChannels.size(); ++i )
                memcpy( outputChannels[i].data( vertIndex ), m_inputChannels[i].data( dataIndex1 ),
                        m_inputChannels[i].primitive_size() );
        } else {
            // tbb::spin_rw_mutex::scoped_lock meshLock(meshMutex, false);
            //  There's no data available, so just add the vertex and set it to zero.
            for( unsigned i = 0; i < outputChannels.size(); ++i )
                memset( outputChannels[i].data( vertIndex ), 0, outputChannels[i].primitive_size() );
        }
    }
}

/**
 *	This function populates all vertex channels specified in the channel propagation
 *	policy using info from the surface policy.
 *	<p>
 *	All other channel info is preserved.
 *
 *	@param	outMesh		Mesh containing vertex and channel info.
 *	@param	cpp			A channel propagation policy.
 */
void direct_linear_rle_level_set_is_policy::populate_mesh_channels(
    frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_propagation_policy& cpp ) const {
    std::vector<frantic::tstring> previousNames; // all channels names in the mesh
    std::vector<frantic::tstring> names;         // all channels that are going to be populated
    std::vector<std::pair<std::size_t, frantic::channels::data_type_t>>
        dataTypes; // types and arity of each channel to be populated
    std::vector<frantic::volumetrics::levelset::const_rle_channel_general_accessor>
        inputChannels; // channel accessors to channel_map
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>
        outputChannels; // channel accessors to vertex_channel

    // If channels already exist throw error.
    outMesh.get_vertex_channel_names( previousNames );
    for( size_t i = 0; i < previousNames.size(); ++i )
        if( cpp.is_channel_included( previousNames[i] ) )
            throw std::runtime_error(
                "direct_linear_rle_level_set_is_policy.populate_mesh_channels: Attempt to populate a "
                "channel that already exists." );

    // Find the names of all channels found in both cpp and the particles.
    get_channel_names( cpp, names );

    // Prepare the input channel accessors and get the data types which
    // will be used to create new channels in the mesh.
    dataTypes.reserve( names.size() );
    inputChannels.reserve( names.size() );
    for( size_t i = 0; i < names.size(); ++i ) {
        inputChannels.push_back( m_levelSet.get_channel_general_accessor( names[i] ) );
        dataTypes.push_back( std::make_pair( inputChannels.back().arity(), inputChannels.back().data_type() ) );
    }

    // Create the new vertex and get accessors to those vertex channels.
    for( size_t i = 0; i < names.size(); ++i ) {
        outMesh.add_vertex_channel_raw( names[i], dataTypes[i].first, dataTypes[i].second );
        outputChannels.push_back( outMesh.get_vertex_channel_general_accessor( names[i] ) );
    }

    // Populate each vertex channel.
    for( size_t i = 0; i < outMesh.vertex_count(); ++i ) {
        populate_vertex_channels( inputChannels, outputChannels, i, outMesh.get_vertex( i ) );
    }
}

/**
 *	This function populates all vertex channels using data from the implicit
 *	surface policy for the input vertex.
 *
 *	@param	inputChannels	Accessors to particle channels used for population.
 *	@param	outputChannels	Accessors to each channel that should be populated.
 *	@param	vert			Location of the vertex.
 *	@param	vertIndex		Index of the vertex.
 */
void direct_linear_rle_level_set_is_policy::populate_vertex_channels(
    const std::vector<frantic::volumetrics::levelset::const_rle_channel_general_accessor>& inputChannels,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert ) const {
    if( inputChannels.size() != outputChannels.size() )
        throw std::runtime_error( "direct_linear_rle_level_set_is_policy.populate_vertex_channels: inputChannels and "
                                  "outputChannels are not the same size" );

    std::vector<float> deltas( 3 );
    std::vector<float> multipliers( 8 );
    std::vector<boost::int32_t> dataIndicies( 8 );

    get_trilerp_indices( m_levelSet.get_rle_index_spec(), m_levelSet.get_voxel_coord_system().get_voxel_coord( vert ),
                         &deltas[0], &dataIndicies[0] );
    get_trilerp_weights( &deltas[0], &multipliers[0] );

    for( size_t i = 0; i < inputChannels.size(); ++i ) {
        inputChannels[i].get_trilinear_interpolated( &multipliers[0], &dataIndicies[0],
                                                     outputChannels[i].data( vertIndex ) );
    }
}

void direct_linear_rle_level_set_is_policy::populate_vertex_channels(
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outputChannels, size_t vertIndex,
    const frantic::graphics::vector3f& vert, vertex_workspace_t& /*workspace*/ ) const {
    populate_vertex_channels( m_inputChannels, outputChannels, vertIndex, vert );
}

} // namespace implicitsurface
} // namespace volumetrics
} // namespace frantic
