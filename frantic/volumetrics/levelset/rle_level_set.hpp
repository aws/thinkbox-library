// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/shared_array.hpp>

#include <frantic/fluids/rle_voxel_field.hpp>
#include <frantic/graphics2d/boundrect2.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec_block_iterator_x.hpp>
#include <frantic/volumetrics/levelset/rle_named_channels.hpp>
#include <frantic/volumetrics/rle_plane.hpp>
#include <frantic/volumetrics/voxel_coord_system.hpp>

#include <frantic/graphics/ray3f.hpp>

#include <frantic/logging/progress_logger.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

/**
 * The RLE Level Set uses an instance of the rle_index_spec to hold a sparsely-stored level set.  It has a primary
 * data array consisting of signed distances to the surface, as well as an arbitrary number of named data channels.
 */
class rle_level_set {
    rle_index_spec m_rleIndex;
    voxel_coord_system m_voxelCoordSystem;
    // TODO: I recommend getting rid of all four of these distances for refactoring rle_level_set into rle_voxel_field.
    //       It seems like just hardcoding 100*voxelLength and -100*voxelLength (or something similar) is perfectly
    //       reasonable, and increases the flexibility a *lot* when the signed distance is then a general float channel
    //       instead of special cased. - Mark W.
    // These are the voxel distance away from the surface towards the inside and the outside.  They are separate,
    // so that we can concentrate the storage on the part (inside or outside) we need more data for.
    float m_interfaceVoxelWidthInside,
        m_interfaceVoxelWidthOutside; // TODO: these are radii or half-widths, rename them?
    // These are the world-space distances to use for undefined voxels that are "inside" or "outside"
    float m_outsideDistance, m_insideDistance;

    // The distance channel, treated specially
    std::vector<float> m_distanceData;

    // Generic named channels
    std::map<frantic::tstring, rle_channel> m_namedChannels;

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
    friend void frantic::fluids::rle_voxel_field::swap_with_rle_level_set( rle_level_set& levelSet );

    friend void detail::build_test1( rle_level_set& levelSet );
    friend void detail::build_testSuite( rle_level_set& levelSet, const frantic::graphics::vector3 orig,
                                         const frantic::graphics::size3 coordSize, const int exteriorRegionCode );
    friend void detail::build_test1_complement( rle_level_set& levelSet );
    friend void detail::build_test2( rle_level_set& levelSet );
    friend void detail::build_test12_union( rle_level_set& levelSet );
    friend void detail::build_test12_intersect( rle_level_set& levelSet );
    friend void detail::build_plane_test( rle_level_set& levelSet );

  public:
    rle_level_set()
        : m_interfaceVoxelWidthInside( 1 )
        , m_interfaceVoxelWidthOutside( 1 )
        , m_outsideDistance( 2 * m_voxelCoordSystem.voxel_length() )
        , m_insideDistance( -2 * m_voxelCoordSystem.voxel_length() ) {}

    rle_level_set( const voxel_coord_system& vcs )
        : m_voxelCoordSystem( vcs )
        , m_interfaceVoxelWidthInside( 1 )
        , m_interfaceVoxelWidthOutside( 1 )
        , m_outsideDistance( 2 * vcs.voxel_length() )
        , m_insideDistance( -2 * vcs.voxel_length() ) {}

    rle_level_set( const voxel_coord_system& vcs, const rle_index_spec& ris, const std::vector<float>& distanceData,
                   float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside )
        : m_rleIndex( ris )
        , m_voxelCoordSystem( vcs )
        , m_interfaceVoxelWidthInside( interfaceVoxelWidthInside )
        , m_interfaceVoxelWidthOutside( interfaceVoxelWidthOutside )
        , m_outsideDistance( ( interfaceVoxelWidthOutside + 1 ) * vcs.voxel_length() )
        , m_insideDistance( -( interfaceVoxelWidthInside + 1 ) * vcs.voxel_length() )
        , m_distanceData( distanceData ) {
        if( distanceData.size() != ris.data_size() )
            throw std::runtime_error( "rle_level_set constructor: The size of the data provided, " +
                                      boost::lexical_cast<std::string>( distanceData.size() ) +
                                      ", doesn't match the data size of the RLE Index Spec provided, " +
                                      boost::lexical_cast<std::string>( ris.data_size() ) + "." );
    }

    /**
     * Returns true if there are not any defined voxels in the level set.
     */
    bool empty() const { return m_distanceData.empty(); }

    // Provide array-based access to the data
    std::size_t size() const { return m_distanceData.size(); }

    void swap( rle_level_set& rhs ) {
        m_rleIndex.swap( rhs.m_rleIndex );
        m_voxelCoordSystem.swap( rhs.m_voxelCoordSystem );
        std::swap( m_interfaceVoxelWidthInside, rhs.m_interfaceVoxelWidthInside );
        std::swap( m_interfaceVoxelWidthOutside, rhs.m_interfaceVoxelWidthOutside );
        std::swap( m_outsideDistance, rhs.m_outsideDistance );
        std::swap( m_insideDistance, rhs.m_insideDistance );
        m_distanceData.swap( rhs.m_distanceData );
        m_namedChannels.swap( rhs.m_namedChannels );
    }

    void clear() {
        m_rleIndex.clear();
        m_distanceData.clear();
        m_namedChannels.clear();
    }

    /**
     * Provide array-based access to the data
     */
    float& operator[]( std::size_t index ) { return m_distanceData[index]; }

    /**
     * Provide array-based access to the data
     */
    const float& operator[]( std::size_t index ) const { return m_distanceData[index]; }

    /**
     * Provide an [] operator using voxel coordinate.  Throws an exception when trying to access invalid data.
     */
    float& operator[]( const frantic::graphics::vector3& voxel ) {
        int dataIndex = m_rleIndex.XYZtoDataIndex( voxel );
        if( dataIndex >= 0 )
            return m_distanceData[dataIndex];
        else
            throw std::runtime_error( "rle_level_set.operator[]: Tried to access undefined voxel coordinate " +
                                      voxel.str() + ", which has region code " +
                                      boost::lexical_cast<std::string>( dataIndex ) );
    }

    /**
     * Provide an [] operator using voxel coordinate.  Throws an exception when trying to access invalid data.
     */
    const float& operator[]( const frantic::graphics::vector3& voxel ) const {
        int dataIndex = m_rleIndex.XYZtoDataIndex( voxel );
        if( dataIndex >= 0 )
            return m_distanceData[dataIndex];
        else
            throw std::runtime_error( "rle_level_set.operator[]: Tried to access undefined voxel coordinate " +
                                      voxel.str() + ", which has region code " +
                                      boost::lexical_cast<std::string>( dataIndex ) );
    }

    /**
     * This retrieves the data from a raw data index.  If the data index is non-negative, it
     * returns the data from the array, otherwise it returns +-m_outsideDistance depending
     * on the region code.
     *
     * @param  dataIndex  The data index to convert to a level set signed distance.
     */
    float get_using_data_index( boost::int32_t dataIndex ) const {
        // casting to avoid having to perform the extra check
        if( (boost::uint32_t)dataIndex < (boost::uint32_t)m_distanceData.size() )
            return m_distanceData[dataIndex];
        else
            return dataIndex == -1 ? m_outsideDistance : m_insideDistance;
    }

    /**
     * This gets the data index to a given voxel, which you can then use with the [] operator to get and set the data.
     * If the voxel is undefined it returns the region code.
     * To test if a voxel is defined you check that the returned value is >= 0.
     *
     * @todo Come up with a good naming scheme!
     *
     * @param  voxel  The voxel coordinate to convert to a data index.
     */
    int XYZtoDataIndex( const frantic::graphics::vector3& voxel ) const { return m_rleIndex.XYZtoDataIndex( voxel ); }

    /**
     * DEPRECATED:
     * This function is to be replaced by frantic::volumetrics::levelset::trilerp_float.
     * This gets the distance from the interface of a given voxel coordinate, using trilinear interpolation.  The value
     * returned is only a valid distance within the interfaceWorldWidth of the interface, assuming the actual data is
     * wide enough and has been appropriately initialized to a distance function.
     */
    float trilerp_voxel_to_signed_distance( const frantic::graphics::vector3f& voxelCoordinate ) const {
        boost::int32_t dataIndices[8];
        // The data samples are at the center of each integer voxel, so we have to offset by 0.5 to compensate.
        float dx = voxelCoordinate.x - 0.5f, dy = voxelCoordinate.y - 0.5f, dz = voxelCoordinate.z - 0.5f;
        frantic::graphics::vector3f voxelMin( floorf( dx ), floorf( dy ), floorf( dz ) );
        dx -= voxelMin.x;
        dy -= voxelMin.y;
        dz -= voxelMin.z;
        // Get all 8 data index values.
        m_rleIndex.fill_2x2x2_data_index_box(
            frantic::graphics::vector3( (int)voxelMin.x, (int)voxelMin.y, (int)voxelMin.z ), dataIndices );

        return ( 1 - dz ) * ( ( 1 - dy ) * ( ( 1 - dx ) * get_using_data_index( dataIndices[0] ) +
                                             dx * get_using_data_index( dataIndices[1] ) ) +
                              dy * ( ( 1 - dx ) * get_using_data_index( dataIndices[2] ) +
                                     dx * get_using_data_index( dataIndices[3] ) ) ) +
               dz * ( ( 1 - dy ) * ( ( 1 - dx ) * get_using_data_index( dataIndices[4] ) +
                                     dx * get_using_data_index( dataIndices[5] ) ) +
                      dy * ( ( 1 - dx ) * get_using_data_index( dataIndices[6] ) +
                             dx * get_using_data_index( dataIndices[7] ) ) );
    }

    /**
     * DEPRECATED:
     * This function is to be replaced by frantic::volumetrics::levelset::trilerp_float.
     */
    float trilerp_signed_distance( const frantic::graphics::vector3f& worldCoordinate ) const {
        return trilerp_voxel_to_signed_distance( m_voxelCoordSystem.get_voxel_coord( worldCoordinate ) );
    }

    /**
     * DEPRECATED:
     * This function is to be replaced by frantic::volumetrics::levelset::trilerp_float.
     * This gets the distance from the interface of a given voxel coordinate, using trilinear interpolation.  The value
     * returned is only a valid distance within the interfaceWorldWidth of the interface, assuming the actual data is
     * wide enough and has been appropriately initialized to a distance function. This function also returns the
     * gradient of the distance field
     */
    float trilerp_voxel_to_signed_distance( const frantic::graphics::vector3f& voxelCoordinate,
                                            frantic::graphics::vector3f& outSignedDistanceGradient ) const {
        boost::int32_t dataIndices[8];
        // The data samples are at the center of each integer voxel, so we have to offset by 0.5 to compensate.
        float dx = voxelCoordinate.x - 0.5f, dy = voxelCoordinate.y - 0.5f, dz = voxelCoordinate.z - 0.5f;
        frantic::graphics::vector3f voxelMin( floorf( dx ), floorf( dy ), floorf( dz ) );
        dx -= voxelMin.x;
        dy -= voxelMin.y;
        dz -= voxelMin.z;
        // Get all 8 data index values.
        m_rleIndex.fill_2x2x2_data_index_box(
            frantic::graphics::vector3( (int)voxelMin.x, (int)voxelMin.y, (int)voxelMin.z ), dataIndices );

        // Get all the samples
        float signedDistance[8];
        for( int i = 0; i < 8; ++i )
            signedDistance[i] = get_using_data_index( dataIndices[i] );

        // Reconstruct the gradient, which is relative to voxel coordinates.  A change of variable to world coordinates
        // requires a scale.
        outSignedDistanceGradient.x = ( 1 - dz ) * ( ( 1 - dy ) * ( signedDistance[1] - signedDistance[0] ) +
                                                     dy * ( signedDistance[3] - signedDistance[2] ) ) +
                                      dz * ( ( 1 - dy ) * ( signedDistance[5] - signedDistance[4] ) +
                                             dy * ( signedDistance[7] - signedDistance[6] ) );
        outSignedDistanceGradient.y = ( 1 - dz ) * ( ( ( 1 - dx ) * signedDistance[2] + dx * signedDistance[3] ) -
                                                     ( ( 1 - dx ) * signedDistance[0] + dx * signedDistance[1] ) ) +
                                      dz * ( ( ( 1 - dx ) * signedDistance[6] + dx * signedDistance[7] ) -
                                             ( ( 1 - dx ) * signedDistance[4] + dx * signedDistance[5] ) );
        outSignedDistanceGradient.z = ( ( 1 - dy ) * ( ( 1 - dx ) * signedDistance[4] + dx * signedDistance[5] ) +
                                        dy * ( ( 1 - dx ) * signedDistance[6] + dx * signedDistance[7] ) ) -
                                      ( ( 1 - dy ) * ( ( 1 - dx ) * signedDistance[0] + dx * signedDistance[1] ) +
                                        dy * ( ( 1 - dx ) * signedDistance[2] + dx * signedDistance[3] ) );

        // Return the trilerped signed distance
        return ( 1 - dz ) * ( ( 1 - dy ) * ( ( 1 - dx ) * signedDistance[0] + dx * signedDistance[1] ) +
                              dy * ( ( 1 - dx ) * signedDistance[2] + dx * signedDistance[3] ) ) +
               dz * ( ( 1 - dy ) * ( ( 1 - dx ) * signedDistance[4] + dx * signedDistance[5] ) +
                      dy * ( ( 1 - dx ) * signedDistance[6] + dx * signedDistance[7] ) );
    }

    /**
     * DEPRECATED:
     * This function is to be replaced by frantic::volumetrics::levelset::trilerp_float.
     * This gets the distance from the interface of a given voxel coordinate, using trilinear interpolation.  The value
     * returned is only a valid distance within the interfaceWorldWidth of the interface, assuming the actual data is
     * wide enough and has been appropriately initialized to a distance function. This function also returns the
     * gradient of the distance field
     */
    float trilerp_signed_distance( const frantic::graphics::vector3f& worldCoordinate,
                                   frantic::graphics::vector3f& outSignedDistanceGradient ) const {
        float result = trilerp_voxel_to_signed_distance( m_voxelCoordSystem.get_voxel_coord( worldCoordinate ),
                                                         outSignedDistanceGradient );
        outSignedDistanceGradient /= m_voxelCoordSystem.voxel_length();
        return result;
    }

    /**
     * Provide an == operator to compare two RLE level sets.
     *
     * @todo  I'm inclined to remove this operator, and have a function similar to this for debugging purposes.  This
     *        shouldn't be a common operation to use.
     */
    bool operator==( const rle_level_set& rls ) {
        if( !( m_voxelCoordSystem == rls.m_voxelCoordSystem ) )
            return false;

        if( !( m_distanceData == rls.m_distanceData ) )
            return false;

        if( !( m_rleIndex == rls.m_rleIndex ) )
            return false;

        //		if( !( m_namedChannels == rls.m_namedChannels ) )
        //			return false;

        return true;
    }

    /**
     * Provide access to the voxel coordinate system.
     */
    const voxel_coord_system& get_voxel_coord_system() const { return m_voxelCoordSystem; }

    const rle_index_spec& get_rle_index_spec() const { return m_rleIndex; }

    frantic::graphics::boundbox3f world_outer_bounds() const {
        return m_voxelCoordSystem.get_world_bounds( m_rleIndex.outer_bounds() );
    }

    float get_interface_voxel_width_inside() const { return m_interfaceVoxelWidthInside; }

    float get_interface_voxel_width_outside() const { return m_interfaceVoxelWidthOutside; }

    float get_outside_distance() const { return m_outsideDistance; }

    float get_inside_distance() const { return m_insideDistance; }

    const std::vector<float>& get_distance_data() const { return m_distanceData; }

    float compute_volume() const;

    ////////////////////////////////////////
    // Named channel methods
    ////////////////////////////////////////

    /**
     * Checks whether this rle level set has the specified named channel.
     *
     * @param  channelName  The name of the channel to check for.
     */
    bool has_channel( const frantic::tstring& channelName ) const {
        return m_namedChannels.find( channelName ) != m_namedChannels.end();
    }

    /**
     * Erases the specified channel from the rle level set.
     *
     * @param  channelName  The name of the channel to erase.
     */
    void erase_channel( const frantic::tstring& channelName ) { m_namedChannels.erase( channelName ); }

    /**
     * This function duplicates the channel named sourceChannelName, copying it into a channel named destChannelName.
     *
     * @param  destChannelName  Where the duplicated channel will end up.
     * @param  sourceChannelName  The source channel to duplicate.
     */
    void duplicate_channel( const frantic::tstring& destChannelName, const frantic::tstring& sourceChannelName );

    /**
     * Copy the level set's signed distance data into a channel named
     * destChannelName.
     *
     * @param  destChannelName  Where the signed distance will be copied to.
     */
    void duplicate_signed_distance_channel( const frantic::tstring& destChannelName );

    /**
     * This retrieves the names of all the channels.  It adds the names to the existing vector, so if you want
     * just the names this function provides, make sure to initialize the output vector to empty yourself.
     *
     * @param  outNames  A vector into which all the names get added.
     */
    void get_channel_names( std::vector<frantic::tstring>& outNames ) const {
        outNames.reserve( outNames.size() + m_namedChannels.size() );
        for( std::map<frantic::tstring, rle_channel>::const_iterator i = m_namedChannels.begin();
             i != m_namedChannels.end(); ++i ) {
            outNames.push_back( i->first );
        }
    }

    /**
     * This adds a channel to the rle level set, of the given name, using the template type parameter for
     * determining the arity and data type.
     *
     * It doesn't initialize the data in the named channel, so if you want all zeros or something like
     * that you have to initialize the data yourself.
     *
     * @param  channelName  The name of the channel to add.
     */
    template <class DataType>
    void add_channel( const frantic::tstring& channelName ) {
        add_channel( channelName, frantic::channels::channel_data_type_traits<DataType>::arity(),
                     frantic::channels::channel_data_type_traits<DataType>::data_type() );
    }

    /**
     * This adds a channel to the rle level set, of the given name, data type and arity.  If this channel already
     * exists, the function does nothing.  If the channel has the wrong data type or arity, it is recreated.
     *
     * It doesn't initialize the data in the named channel, so if you want all zeros or something like
     * that you have to initialize the data yourself.
     *
     * @param  channelName  The name of the channel to add.
     * @param  arity  The arity of the channel.
     * @param  dataType  The data type of the channel.
     */
    void add_channel( const frantic::tstring& channelName, std::size_t arity, data_type_t dataType );

    /**
     * This sets the specified channel to all zeros.
     *
     * @param  channelName  The name of the channel to set to all zeros.
     */
    void zero_channel( const frantic::tstring& channelName );

    /**
     * This function gets the maximum L2-norm of a value in the given channel.
     *
     * @param  channelName  The name of the channel to set to all zeros.
     */
    float get_channel_max_norm( const frantic::tstring& channelName ) const;

    /**
     * This function copies the specified channels from the input level set.  Any defined voxels in the destination
     * which are not in the input are set to have channel value 0.  If requested, an uint8 named channel called
     * populatedChannelToCreate is created, which indicates where values were set.  When a 'Populated' channel already
     * exists, and is requested by the caller, the data within the existing channel gets overwritten.
     *
     * @param  inputRLS  The level set from which to copy the channels.
     * @param  channelsToCopy  A list of channel names which should be copied.
     * @param  populatedChannelToCreate  If non-empty, the function adds a new uint8 named channel with the given name,
     * as a 'Populated' channel.
     */
    void copy_channels( const rle_level_set& inputRLS, const std::vector<frantic::tstring>& channelsToCopy,
                        const frantic::tstring& populatedChannelToCreate = _T("") );

    void copy_channels( const frantic::fluids::rle_voxel_field& inputRVF,
                        const std::vector<frantic::tstring>& channelsToCopy,
                        const frantic::tstring& populatedChannelToCreate = _T("") );

    rle_channel_general_accessor get_channel_general_accessor( const frantic::tstring& channelName ) {
        std::map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.find( channelName );
        if( i != m_namedChannels.end() ) {
            return i->second.get_general_accessor();
        } else {
            throw std::runtime_error(
                "rle_level_set.get_channel_general_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( channelName ) + "\", but no such channel exists in this rle_level_set." );
        }
    }

    const_rle_channel_general_accessor get_channel_general_accessor( const frantic::tstring& channelName ) const {
        std::map<frantic::tstring, rle_channel>::const_iterator i = m_namedChannels.find( channelName );
        if( i != m_namedChannels.end() ) {
            return i->second.get_general_accessor();
        } else {
            throw std::runtime_error(
                "rle_level_set.get_channel_general_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( channelName ) + "\", but no such channel exists in this rle_level_set." );
        }
    }

    template <class DataType>
    rle_channel_accessor<DataType> get_channel_accessor( const frantic::tstring& channelName ) {
        std::map<frantic::tstring, rle_channel>::iterator i = m_namedChannels.find( channelName );
        if( i != m_namedChannels.end() ) {
            return i->second.get_accessor<DataType>();
        } else {
            throw std::runtime_error(
                "rle_level_set.get_channel_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( channelName ) + "\", but no such channel exists in this rle_level_set." );
        }
    }

    template <class DataType>
    const_rle_channel_accessor<DataType> get_channel_accessor( const frantic::tstring& channelName ) const {
        std::map<frantic::tstring, rle_channel>::const_iterator i = m_namedChannels.find( channelName );
        if( i != m_namedChannels.end() ) {
            return i->second.get_accessor<DataType>();
        } else {
            throw std::runtime_error(
                "rle_level_set.get_channel_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( channelName ) + "\", but no such channel exists in this rle_level_set." );
        }
    }

    ////////////////////////////////////////
    // Queries to get data
    ////////////////////////////////////////

    /**
     * This function fills a dense xy plane as specified by the bounding rectangle with sample values from the level
     * set.
     *
     * Before calling this method, outVoxelCornerValues must be allocated
     * to hold at least voxelXYExtents.get_area() floats.
     */
    void fill_plane( const frantic::graphics2d::boundrect2& voxelXYExtents, int voxelZ,
                     float* outVoxelCornerValues ) const;

    /**
     * This function fills a dense xy plane as specified by the bounding rectangle with sample values from the level
     * set.
     */
    void fill_plane( const frantic::graphics2d::boundrect2& voxelXYExtents, int voxelZ,
                     std::vector<float>& outVoxelCornerValues ) const;

    /**
     * This function fills a sparse xy plane as specified by the bounding rectangle with sample values from the level
     * set.
     */
    //	void fill_sparse_plane( const frantic::graphics2d::boundrect2& voxelXYExtents,
    //							int voxelZ,
    //							float* outVoxelCornerValues,
    //							frantic::volumetrics::rle_plane& outRLP ) const;

    /**
     * This function fills sparse xy planes as specified by the bounding rectangle with data from a level set
     */
    void fill_sparse_plane_channel_data( const frantic::graphics2d::boundrect2& xyExtents, int z,
                                         std::vector<frantic::tstring>& channelNames, std::vector<char*>& channelData,
                                         frantic::volumetrics::rle_plane& rlp ) const;

    /**
     * This function fills a dense integer box with sample values from the level set.
     *
     * Before calling this method, outVoxelCornerValues must be allocated
     * to hold at least voxelExtents.get_volume() floats.
     */
    void fill_box( const frantic::graphics::boundbox3& voxelExtents, float* outVoxelCornerValues ) const;

    /**
     * This function fills a dense integer box with sample values from the level set.
     */
    void fill_box( const frantic::graphics::boundbox3& voxelExtents, std::vector<float>& outVoxelCornerValues ) const;

    /**
     * This function fills a dense xy plane as specified by the bounding rectangle with sample values from the level
     * set. It uses a reconstruction filter to provide the values on the provided coordinate system.
     *
     * @todo  Add selection of the reconstruction filter as well.
     */
    template <class ReconstructionFilter>
    void fill_plane( const voxel_coord_system& destCoordSys, const ReconstructionFilter& reconFilter,
                     const frantic::graphics2d::boundrect2& voxelXYExtents, int voxelZ,
                     std::vector<float>& outVoxelCornerValues ) const;
    template <class ReconstructionFilter>
    void fill_plane_mt( const voxel_coord_system& destCoordSys, const ReconstructionFilter& reconFilter,
                        const frantic::graphics2d::boundrect2& voxelXYExtents, int voxelZ,
                        float* outVoxelCornerValues ) const;

    /**
     *	This function fills a z plane with sparse data as indicated by the defined regions of the level set.
     *	It also populates an accompanying rle plane structure for accessing the sparse data.
     *
     *	@param	destCoordSys	The voxel coord system that the plane to be populated lives in.  Doesn't have to
     *							be the same as the level set coord system.  The templated filter will be
     *used to resample the data.
     *	@param	reconFilter		The ilter to be used to reconstruct data values (eg math::bspline).
     *	@param	voxelExtents	The extents of the plane to be populated in the xy direction
     *	@param	voxelZ			The z coord of the plane.
     *	@param	channelNames	Names of channels to fetch.
     *	@param	channelData		This is where the channel data will get put.
     *	@param	outRLP			An rle plane structure for accessing the data sparsely.
     *
     *  @todo  Add selection of the reconstruction filter as well.
     */
    template <class ReconstructionFilter>
    void fill_sparse_plane_channel_data( const voxel_coord_system& destCoordSys,
                                         const ReconstructionFilter& reconFilter,
                                         const frantic::graphics2d::boundrect2& voxelExtents, int voxelZ,
                                         std::vector<frantic::tstring>& channelNames, std::vector<char*>& channelData,
                                         frantic::volumetrics::rle_plane& outRLP ) const;

    template <class ReconstructionFilter>
    void fill_sparse_plane_channel_data( const voxel_coord_system& destCoordSys,
                                         const ReconstructionFilter& reconFilter,
                                         const frantic::graphics2d::boundrect2& voxelExtents, int voxelZ,
                                         std::vector<frantic::tstring>& channelNames,
                                         std::vector<char*>& channelData ) const;

    /**
     * This function fills a dense integer box with sample values from the level set.  It uses a reconstruction filter
     * to provide the values on the provided coordinate system.
     *
     * @todo  Add selection of the reconstruction filter as well.
     */
    template <class ReconstructionFilter>
    void fill_box( const voxel_coord_system& destCoordSys, const ReconstructionFilter& reconFilter,
                   const frantic::graphics::boundbox3& voxelExtents, std::vector<float>& outVoxelCornerValues ) const;
    template <class ReconstructionFilter>
    void fill_box_mt( const voxel_coord_system& destCoordSys, const ReconstructionFilter& reconFilter,
                      const frantic::graphics::boundbox3& voxelExtents, float* outVoxelCornerValues ) const;

    // If all the voxels touching the box have the same region code, return it, otherwise return 0.
    // This function can be used to tell whether a specified box region can be skipped because it's entirely outside,
    // processed trivially because it's entirely inside, or needs to be subdivided because there's likely defined
    // data inside.
    int get_unique_region_code( const frantic::graphics::boundbox3f& worldBounds ) const {
        return m_rleIndex.get_unique_region_code( m_voxelCoordSystem.get_voxel_bounds( worldBounds ) );
    }

    // This prints out an rle level set in 2D
    void print2d( std::ostream& out ) const;

    /*
     * Dumps the data of the rle_level_set to the provided output stream.
     *
     * @param   out  The ostream the structure is dumped to.
     */
    void dump( std::ostream& out ) const;

    ////////////////////////////////////////
    // Functions which construct an rle level set
    ////////////////////////////////////////

    /**
     * This is a convenience function designed to make it easy to create a level set.
     *
     * @note Do not use this function in performance critical code.
     *
     * @param  values  This std::map maps vector3 voxel coordinates into level set distance values.
     * @param  interfaceVoxelWidthInside  This is the value to assume for the inside defined voxel width.
     * @param  interfaceVoxelWidthOutside  This is the value to assume for the outside defined voxel width.
     */
    void set( const std::map<frantic::graphics::vector3, float>& values, float interfaceVoxelWidthInside,
              float interfaceVoxelWidthOutside );

    /**
     * This sets the rle level set to the given data.
     *
     * @param  vcs  The voxel coordinate system to use.
     * @param  ris  The rle index spec to use.
     * @param  distanceData  An array of level set signed distances.  Its size must match ris.data_size().
     * @param  interfaceVoxelWidthInside  This is the value to assume for the inside defined voxel width.
     * @param  interfaceVoxelWidthOutside  This is the value to assume for the outside defined voxel width.
     */
    void set( const voxel_coord_system& vcs, const rle_index_spec& ris, const std::vector<float>& distanceData,
              float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside );

    /**
     * This sets the rle level set to an empty run, defaulting to an "outside" exterior region code.
     *
     * @param  exteriorRegionCode  The region code to use for voxels that don't have a value specified by the
     * rle_index_spec.
     */
    void set_to_empty( boost::int32_t exteriorRegionCode = -1 );

    /**
     * This sets the rle level set to the given data, using swap functions to set the data efficiently.  Note that
     * the variables you pass in will have arbitrary values in them after this function is called.
     *
     * @param  vcs  The voxel coordinate system to use.
     * @param  ris  The rle index spec to use.  Note that the parameter you pass in will be destroyed, because it is
     * swapped into the rle_level_set to avoid the performance and memory overhead of doing a copy.
     * @param  distanceData  An array of level set signed distances.  Its size must match ris.data_size().
     * @param  interfaceVoxelWidthInside  This is the value to assume for the inside defined voxel width.
     * @param  interfaceVoxelWidthOutside  This is the value to assume for the outside defined voxel width.
     */
    void set_with_swap( voxel_coord_system& vcs, rle_index_spec& ris, std::vector<float>& distanceData,
                        float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside );

    /**
     * This function switches the rle_index_spec being used by this rle_level_set.  It sets any new named channel values
     * to 0, and previously undefined level set distance values to m_outsideDistance or m_insideDistance.  If you ask it
     * to, it will also create a named channel called 'Populated' which is a uint8 channel that has 1 for values that
     * already had data and 0 for values that got new data.  This channel can then later be used for a partial
     * fast-marching method processing step which solves for the level set distance values and extrapolates named
     * channel values.
     *
     * @param  ris  The rle_index_spec to switch to.
     * @param  populatedChannelToCreate  If non-empty, the function adds a new uint8 named channel with the given name,
     * as a 'Populated' channel.
     */
    void switch_rle_index_spec( const rle_index_spec& ris, const frantic::tstring& populatedChannelToCreate = _T("") ) {
        rle_index_spec rleTemp = ris;
        switch_rle_index_spec_with_swap( rleTemp, populatedChannelToCreate );
    }

    /**
     * This function switches the rle_index_spec being used by this rle_level_set.  It sets any new named channel values
     * to 0, and previously undefined level set distance values to m_outsideDistance or m_insideDistance.  If you ask it
     * to, it will also create a named channel called 'Populated' which is a uint8 channel that has 1 for values that
     * already had data and 0 for values that got new data.  This channel can then later be used for a partial
     * fast-marching method processing step which solves for the level set distance values and extrapolates named
     * channel values.
     *
     * Use this function if you don't need the data in the rle_index_spec you'	re switching it to, and you will avoid
     * some memory allocation and copying.
     *
     * @param  ris  The rle_index_spec to switch to.  The data of this ris is stolen and used in the rle_level_set, to
     *              minimize memory overhead.
     * @param  populatedChannelToCreate  If non-empty, the function adds a new uint8 named channel with the given name,
     * as a 'Populated' channel.
     */
    void switch_rle_index_spec_with_swap( rle_index_spec& ris,
                                          const frantic::tstring& populatedChannelToCreate = _T("") );

    /**
     * Creates an int32 channel which maps each defined voxel in this level set to the data index of the defined voxel
     * for the destination rle_index_spec.
     *
     * @param  dataIndexMapChannelToCreate  The name of the channel to create
     * @param  risMappingTarget  The rle_index_spec which is the target to which the created channel refers.
     */
    void create_data_index_map_channel( const frantic::tstring& dataIndexMapChannelToCreate,
                                        const rle_index_spec& risMappingTarget );

    /**
     * This function creates a centered volume channel and/or a staggered volume channel, using this level set for the
     * input.
     *
     * The centered volume channel consists of the voxel-normalized volume of the fluid within each voxel.  Its value
     * ranges from 0.0 to 1.0.
     *
     * What the staggered volume channel means, is that for each face, the field has the volume within the voxel-sized
     * cube centered about that face, relative to the voxel size. This means for entirely inside regions, the values
     * will be 1.0, and for entirely outside regions, the values will be 0.0.
     *
     * The centered channel created has type float32, and the staggered channel created has type float32[3].
     *
     * @param  outputField  The rle_voxel_field within which the volume channel(s) will be created.
     * @param  outputCenteredVolumeChannelName  The name of the centered volume channel to be created.  Specify "" to
     * skip this channel.
     * @param  outputStaggeredVolumeChannelName  The name of the staggered volume channel to be created.  Specify "" to
     * skip this channel.
     */
    void create_volume_fraction_channels( fluids::rle_voxel_field& outputField,
                                          const frantic::tstring& outputCenteredVolumeChannelName,
                                          const frantic::tstring& outputStaggeredVolumeChannelName ) const;

    /**
     * This function creates a centered volume channel and/or a staggered volume channel, using this level set for the
     * input.
     *
     * The centered volume channel consists of the voxel-normalized volume of the fluid within each voxel.  Its value
     * ranges from 0.0 to 1.0.
     *
     * What the staggered volume channel means, is that for each face, the field has the volume within the voxel-sized
     * cube centered about that face, relative to the voxel size. This means for entirely inside regions, the values
     * will be 1.0, and for entirely outside regions, the values will be 0.0.
     *
     * The centered channel created has type float32, and the staggered channel created has type float32[3].
     *
     * @param  outputLS  The rle_level_set within which the volume channel(s) will be created.
     * @param  outputCenteredVolumeChannelName  The name of the centered volume channel to be created.  Specify "" to
     * skip this channel.
     * @param  outputStaggeredVolumeChannelName  The name of the staggered volume channel to be created.  Specify "" to
     * skip this channel.
     */
    void create_volume_fraction_channels( rle_level_set& outputLS,
                                          const frantic::tstring& outputCenteredVolumeChannelName,
                                          const frantic::tstring& outputStaggeredVolumeChannelName ) const;

    /**
     * Advects the distance data of the level set using the staggered velocity channel of the voxel field. This is a
     * singled threaded implementation
     *
     * @param  velocityField  The voxel field that contains the Staggered Velocities to use
     * @param  staggeredVelocityChannelName  The name of the staggered velocity channel.
     * @param  dt  The time step for advection.
     */
    void serial_semi_lagrangian_advect_staggered( const frantic::fluids::rle_voxel_field& velocityField,
                                                  const frantic::tstring& staggeredVelocityChannelName, float dt );

    /**
     * Advects the distance data of the level set using the staggered velocity channel of the voxel field.
     * This overload uses the voxelBounds parameter to handle half-open liquid surface level sets.
     *
     * @param  voxelBounds  The simulation bounds, which are used to handle half-open level sets.
     * @param  velocityField  The voxel field that contains the Staggered Velocities to use
     * @param  staggeredVelocityChannelName  The name of the staggered velocity channel.
     * @param  dt  The time step for advection.
     */
    void semi_lagrangian_advect_staggered( const frantic::graphics::boundbox3& voxelBounds,
                                           const frantic::fluids::rle_voxel_field& velocityField,
                                           const frantic::tstring& staggeredVelocityChannelName, float dt );

    /**
     * Advects the distance data of the level set using the staggered velocity channel of the voxel field.
     *
     * @param  velocityField  The voxel field that contains the Staggered Velocities to use
     * @param  staggeredVelocityChannelName  The name of the staggered velocity channel.
     * @param  dt  The time step for advection.
     */
    void semi_lagrangian_advect_staggered( const frantic::fluids::rle_voxel_field& velocityField,
                                           const frantic::tstring& staggeredVelocityChannelName, float dt );

    /**
     * Advects the distance data of the level set using the velocity channel of the voxel field.
     *
     * @param  velocityField  The voxel field that contains the Staggered Velocities to use
     * @param  velocityChannelName  The name of the velocity channel.
     * @param  dt  The time step for advection.
     */
    void semi_lagrangian_advect( const rle_level_set& velocityField, const frantic::tstring& velocityChannelName,
                                 float dt );

    /**
     * This function copies the specified channels from the input level set, while doing a semi-lagrangian advection.
     * It uses the staggered velocity field of the provided voxel field.  If a 'Populated' channel is requested, it is
     * set to 1 whenever the semi-lagrangian query found a fully defined value, and 0 when any of the voxels that went
     * into the interpolated value were undefined.
     *
     * @param  inputRLS  The level set from which to advect and copy the channels.
     * @param  inputStaggeredVelocityField the field that contains a Staggered Velocity channel to use for the advection
     * @param  timeStep  This is the length of the step (in seconds, because the Velocity channel will be world
     * units/second) for the advection.
     * @param  channelsToAdvect  A list of channel names which should be advected and copied.
     * @param  targetVoxelsChannelName  If non-empty, the function expects this to be a uint8 channel with 0 indicating
     * not to target this voxel with a semi-lagrangian advection
     * @param  populatedChannelName  If non-empty, the function adds a new uint8 named channel with 'Populated' data in
     * this channel.
     */
    void semi_lagrangian_advect_channels_staggered( const rle_level_set& inputRLS,
                                                    const fluids::rle_voxel_field& inputStaggeredVelocityField,
                                                    float timeStep,
                                                    const std::vector<frantic::tstring>& channelsToAdvect,
                                                    const frantic::tstring& targetVoxelsChannelName = _T(""),
                                                    const frantic::tstring& populatedChannelName = _T("") );

    /**
     * This function copies the specified channels from the input level set, while doing a semi-lagrangian advection.
     * It uses the staggered velocity field of the provided voxel field.  If a 'Populated' channel is requested, it is
     * set to 1 whenever the semi-lagrangian query found a fully defined value, and 0 when any of the voxels that went
     * into the interpolated value were undefined.
     *
     * @param  inputRLS  The level set from which to advect and copy the channels.
     * @param  inputStaggeredVelocityField the field that contains a Staggered Velocity channel to use for the advection
     * @param  timeStep  This is the length of the step (in seconds, because the Velocity channel will be world
     * units/second) for the advection.
     * @param  channelsToAdvect  A list of channel names which should be advected and copied.
     * @param  targetVoxelsChannelName  If non-empty, the function expects this to be a uint8 channel with 0 indicating
     * not to target this voxel with a semi-lagrangian advection
     * @param  populatedChannelName  If non-empty, the function adds a new uint8 named channel with 'Populated' data in
     * this channel.
     */
    void serial_semi_lagrangian_advect_channels_staggered( const rle_level_set& inputRLS,
                                                           const fluids::rle_voxel_field& inputStaggeredVelocityField,
                                                           float timeStep,
                                                           const std::vector<frantic::tstring>& channelsToAdvect,
                                                           const frantic::tstring& targetVoxelsChannelName = _T(""),
                                                           const frantic::tstring& populatedChannelName = _T("") );

    /**
     * This function copies the specified channels from the input level set, while doing a semi-lagrangian advection.
     * It uses the velocity field of the target level set.  If a 'Populated' channel is requested, it is set to 1
     * whenever the semi-lagrangian query found a fully defined value, and 0 when any of the voxels that went into
     * the interpolated value were undefined.
     *
     * @param  inputRLS  The level set from which to advect and copy the channels.
     * @param  timeStep  This is the length of the step (in seconds, because the Velocity channel will be world
     * units/second) for the advection.
     * @param  channelsToAdvect  A list of channel names which should be advected and copied.
     * @param  targetVoxelsChannelName  If non-empty, the function expects this to be a uint8 channel with 0 indicating
     * not to update a voxel with a semi-lagrangian advection, or a 1 indicating to update a voxel.
     * @param  populatedChannelName  If non-empty, the function adds a new uint8 named channel with 'Populated' data in
     * this channel.
     */
    void semi_lagrangian_advect_channels( const rle_level_set& inputRLS, float timeStep,
                                          const std::vector<frantic::tstring>& channelsToAdvect,
                                          const frantic::tstring& targetVoxelsChannelName = _T(""),
                                          const frantic::tstring& populatedChannelName = _T("") );

    /**
     * This function copies the specified channels from the input level set, much like semi_lagrangian_advect_channels,
     * however, when the semi-lagrangian query has an undefined value in any of the voxels, we simply qeury the
     * non-advected point in the default data level set. If the source level set has any undefined values, the
     * semi-lagrangian query is done the same as in semi_lagrangian_advect_channels.
     *
     * @param  inputRLS  The level set from which to advect and copy the channels.
     * @param  defaultRLS  The level set from which to copy the channels if the inputRLS is undefined.
     * @param  timeStep  This is the length of the step (in seconds, because the Velocity channel will be world
     * units/second) for the advection.
     * @param  channelsToAdvect  A list of channel names which should be advected and copied.
     * @param  targetVoxelsChannelName  If non-empty, the function expects this to be a uint8 channel with 0 indicating
     * not to update a voxel with a semi-lagrangian advection, or a 1 indicating to update a voxel.
     * @param  populatedChannelName  If non-empty, the function adds a new uint8 named channel with 'Populated' data in
     * this channel.
     */
    void semi_lagrangian_advect_channels_with_undefined_default(
        const rle_level_set& inputRLS, const rle_level_set& defaultRLS, float timeStep,
        const std::vector<frantic::tstring>& channelsToAdvect, const frantic::tstring& targetVoxelsChannelName = _T(""),
        const frantic::tstring& populatedChannelName = _T("") );

    /**
     * This function uses marching extrapolation to extrapolate the specified named channels.  It
     * uses and modifies the 'Populated' channel to do this.  If there is no 'Populated' channel,
     * it throws an exception.
     *
     * To specify that a voxel has correct values, and is thus an extrapolation source, set its 'Populated'
     * channel value to 1.  To specify that a voxel has invalid values, and should be an extrapolation target,
     * set its 'Populated' channel value to 0.  To specify that a voxel has invalid values, but shouldn't
     * be used as an extrapolation target, set its 'Populated' channel value to 2.  All other 'Populated' channel
     * values are reserved.
     *
     * On exit from this function, some voxels may have 'Populated' channel values of 0 or 3.  If the value is 0,
     * that means that the voxel is completely disconnected from any marching sources.  If the value is 3, that
     * means that either there were no valid inputs to extrapolate to this voxel, or the inputs were ambiguous.
     * One ambiguous case is extrapolating across the level set interface, where there are unpopulated voxels
     * on both sides of the interface, creating an extrapolation dependency loop.
     *
     * @note This function allocates and uses the cached ris_adjacency in the rle_index_spec.
     *
     * @param  channelsToExtrapolate  A list of channel names which should be extrapolated.
     * @param  correspondingGradientChannelsToMatch  An list of the same size as channelsToExtrapolate, which specifies
     * gradient channel names to match for the channels being extrapolated.  If an empty string is specified, the
     * default constant extrapolation (equivalent to matching a 0 gradient) is used for the corresponding extrapolation
     * channel.
     * @param  populatedChannelName  The name of the 'Populated' channel, lets you customize/use multiple such channels
     * at once.
     */
    void extrapolate_channels( const std::vector<frantic::tstring>& channelsToExtrapolate,
                               const std::vector<frantic::tstring>& correspondingGradientChannelsToMatch,
                               const frantic::tstring& populatedChannelName );

    /**
     *  Extrapolate the specified named channel, while matching a specified gradient channel.
     *
     * @see void extrapolate_channels( const std::vector<std::string>& channelsToExtrapolate, const
     * std::vector<std::string>& correspondingGradientChannelsToMatch, const std::string& populatedChannelName )
     *
     * @param channelToExtrapolate The name of the channel which should be extrapolated.
     * @param gradientChannelToMatch The name of the gradient channel to match when extrapolating.
     * @param populatedChannelName The name of the 'Populated' channel, lets you customize/use multiple such channels at
     * once.
     */
    void extrapolate_channel( const frantic::tstring& channelToExtrapolate,
                              const frantic::tstring& gradientChannelToMatch,
                              const frantic::tstring& populatedChannelName );

    /**
     * This function uses marching extrapolation to extrapolate the specified named channels.  It
     * uses and modifies the 'Populated' channel to do this.  If there is no 'Populated' channel,
     * it throws an exception.
     *
     * To specify that a voxel has correct values, and is thus an extrapolation source, set its 'Populated'
     * channel value to 1.  To specify that a voxel has invalid values, and should be an extrapolation target,
     * set its 'Populated' channel value to 0.  To specify that a voxel has invalid values, but shouldn't
     * be used as an extrapolation target, set its 'Populated' channel value to 2.  All other 'Populated' channel
     * values are reserved.
     *
     * On exit from this function, populated voxels will have 'Populated' channel value 0, while unpopulated
     * voxels will have 'Populated' channel values of 0 or 3.  If the value is 0,
     * that means that the voxel is completely disconnected from any marching sources.  If the value is 3, that
     * means that either there were no valid inputs to extrapolate to this voxel, or the inputs were ambiguous.
     * One ambiguous case is extrapolating across the level set interface, where there are unpopulated voxels
     * on both sides of the interface, creating an extrapolation dependency loop.
     *
     * @note This function allocates and uses the cached ris_adjacency in the rle_index_spec.
     *
     * @param  channelsToExtrapolate  A list of channel names which should be extrapolated.
     * @param  populatedChannelName  The name of the 'Populated' channel, lets you customize/use multiple such channels
     * at once.
     */
    void extrapolate_channels( const std::vector<frantic::tstring>& channelsToExtrapolate,
                               const frantic::tstring& populatedChannelName );

    /**
     *  Extrapolate the specified named channel.
     *
     * @see void extrapolate_channels( const std::vector<std::string>& channelsToExtrapolate, const std::string&
     * populatedChannelName )
     *
     * @param channelToExtrapolate The name of the channel which should be extrapolated.
     * @param populatedChannelName The name of the 'Populated' channel, lets you customize/use multiple such channels at
     * once.
     */
    void extrapolate_channel( const frantic::tstring& channelToExtrapolate,
                              const frantic::tstring& populatedChannelName );

    /**
     * This extrapolates all the named channels in the level set except the 'Populated' channel.  See the other
     * extrapolate_channels function for more details.
     *
     * @note This function allocates and uses the cached ris_adjacency in the rle_index_spec.
     *
     * @param  populatedChannelName  The name of the 'Populated' channel, lets you customize/use multiple such channels
     * at once.
     */
    void extrapolate_channels( const frantic::tstring& populatedChannelName );

    /**
     * This extrapolates a staggered channel in a separate rle_voxel_field, using the distance data
     * from the rle_level_set self to guide the extrapolation.
     *
     * This does not change the populated values stored in rlsStaggeredPopulatedChannelName.
     *
     * @param  rvf                        The rle_voxel_field which contains the staggered field to extrapolate.
     * @param  rvfStaggeredChannelName    The name of the staggered field channel in the rvf.
     * @param  extrapDirection            If +1, extrapolate from inside to outside, -1, extrapolate from outside to
     * inside.
     * @param  rlsStaggeredPopulatedChannelName  The name of the staggered populated channel in the rls.  Its type
     * should be uint8, and bit 0 is for X, 1 for Y, 2 for Z.
     * @param  rlsDataIndexMapChannelName  If non-empty, this should reference a channel created with the
     *                                     create_data_index_map_channel function.
     */
    void extrapolate_staggered_channel( fluids::rle_voxel_field& rvf, const frantic::tstring& rvfStaggeredChannelName,
                                        int extrapDirection, const frantic::tstring& rlsStaggeredPopulatedChannelName,
                                        const frantic::tstring& rlsDataIndexMapChannelName = _T("") );

    /**
     * This function sets the 'Populated' channel to have a value 1 for all boundary voxels, and a value 0 for
     * all non-boundary voxels.  A boundary voxel is defined as one where the level set signed distance changes
     * sign from one of its neighbor voxels.
     *
     * @note This function allocates and uses the cached ris_adjacency in the rle_index_spec.
     *
     * @param  populatedChannelName  The name of the 'Populated' channel, lets you customize/use multiple such channels
     * at once.
     */
    void tag_interface_voxels( const frantic::tstring& populatedChannelName );

    /**
     * This function sets the 'Populated' channel to have a value 1 for all voxels with a distance value between
     * minVoxelDistance and maxVoxelDistance, and 0 for all other voxels
     *
     * @param  populatedChannelName  The name of the 'Populated' channel, lets you customize/use multiple such channels
     * at once.
     * @param  minVoxelDistance  A signed distance value, in voxel units, that is the minimum distance value to tag.
     * @param  maxVoxelDistance  A signed distance value, in voxel units, that is the maximum distance value to tag.
     */
    void tag_distance_range( const frantic::tstring& populatedChannelName, float minVoxelDistance,
                             float maxVoxelDistance );

    /**
     * This function sets the 'StaggeredPopulated' channel to have a value 1 for all voxel faces with a distance value
     * between minDistance and maxDistance, and 0 for all other voxels.  For the X face, it uses bit mask 0x01, for the
     * Y face, it uses bit max 0x02, and for the Z face, it uses bit mask 0x04.
     *
     * @param  staggeredPopulatedChannelName  The name of the 'StaggeredPopulated' channel, lets you customize/use
     * multiple such channels at once.
     * @param  minVoxelDistance  A signed distance value, in voxel units, that is the minimum distance value to tag.
     * @param  maxVoxelDistance  A signed distance value, in voxel units, that is the maximum distance value to tag.
     */
    void tag_distance_range_staggered( const frantic::tstring& staggeredPopulatedChannelName, float minVoxelDistance,
                                       float maxVoxelDistance );

    /**
     * This function uses the fast marching method to reinitialize the signed distance channels.  It starts marching from
     * the existing isosurface, first updating the values on that band before marching outwards.  If
     * requested, a 'Populated' channel can be created.  On return, this 'Populated' channel will
     * have 0 for unpopulated voxels, and 1 for populated voxels that were successfully reinitialized.
     *
     * @note This function allocates and uses the cached ris_adjacency in the rle_index_spec.
     *
     * @param  logger		the logger will be used to display progress
     * @param  populatedChannelToCreate  The name of the 'Populated' channel that will get created.  By default none is
     * created.
     * @param  insideVoxelDistance  The signed distance value, in voxel units, at which to stop marching on the inside
     * defined region of the level set. This value must be negative.
     * @param  outsideVoxelDistance  The signed distance value, in voxel units, at which to stop marching on the outside
     * defined region of the level set. This value must be positive.
     */
    void reinitialize_signed_distance( frantic::logging::progress_logger& logger,
                                       const frantic::tstring& populatedChannelToCreate,
                                       const float insideVoxelDistance, const float outsideVoxelDistance );

    /// see above
    inline void
    reinitialize_signed_distance( const frantic::tstring& populatedChannelToCreate = _T(""),
                                  const float insideVoxelDistance = -std::numeric_limits<float>::infinity(),
                                  const float outsideVoxelDistance = std::numeric_limits<float>::infinity() ) {

        frantic::logging::null_progress_logger nullProgressLogger;
        reinitialize_signed_distance( nullProgressLogger, populatedChannelToCreate, insideVoxelDistance,
                                      outsideVoxelDistance );
    }

    /**
     * This function uses the fast marching method to reinitialize the signed distance channels.
     * It starts marching from
     * the boundary of the provided 'Populated' channel, and modifies the 'Populated' channel to include
     * the newly filled voxels.  If there is no 'Populated' channel, it throws an exception.
     *
     * To specify that a voxel has correct values, and is thus a reinitialization source, set its 'Populated'
     * channel value to 1.  To specify that a voxel has invalid values, and should be an reinitialization target,
     * set its 'Populated' channel value to 0.  To specify that a voxel has invalid values, but shouldn't
     * be used as a reinitialization target, set its 'Populated' channel value to 2.  All other 'Populated' channel
     * values are reserved.  On return, the 'Populated' channel will be 0 for
     * unpopulated voxels, and 1 for voxels that were successfully reinitialized
     * or were initially populated.
     *
     * @note This function allocates and uses the cached ris_adjacency in the rle_index_spec.
     *
     * @param  progressLogger	 the logger will be used to display progress
     * @param  populatedChannelName  The name of the 'Populated' channel, from which marching starts.
     * @param  insideVoxelDistance  The signed distance value at which to stop marching on the inside defined region of
     * the level set. This value must be negative.
     * @param  outsideVoxelDistance  The signed distance at which to stop marching on the outside defined region of the
     * level set. This value must be positive.
     */
    void reinitialize_signed_distance_from_populated( frantic::logging::progress_logger& progressLogger,
                                                      const frantic::tstring& populatedChannelName,
                                                      const float insideVoxelDistance,
                                                      const float outsideVoxelDistance );

    /// see above
    inline void reinitialize_signed_distance_from_populated(
        const frantic::tstring& populatedChannelName,
        const float insideVoxelDistance = -std::numeric_limits<float>::infinity(),
        const float outsideVoxelDistance = std::numeric_limits<float>::infinity() ) {

        frantic::logging::null_progress_logger nullProgressLogger;
        reinitialize_signed_distance_from_populated( nullProgressLogger, populatedChannelName, insideVoxelDistance,
                                                     outsideVoxelDistance );
    }

  private:
    /**
     * This function implements a CSG operation based on the provided policy class.  It is used for implementing
     * union, intersection, and subtraction.
     *
     * This is an internal, private function.
     *
     * @param   rleFirst                The first rle_level_set operand.
     * @param   rleSecond               The second rle_level_set operand.
     * @param   channelBlendVoxelDistance    The distance, specified in units of voxels, over which the named channel
     * values are blended. A value of 1 or so is usually reasonable.
     */
    template <class CSGPolicy>
    void csg_operation( const rle_level_set& rleFirst, const rle_level_set& rleSecond,
                        float channelBlendVoxelDistance );

  public:
    /**
     * This function converts an rle level set to its complement.  The operation is done in place.
     */
    void csg_complement();

    /**
     * This function overwrites the current rle level set with the CSG union of the two operands.  This is done by
     * taking the minimum value of the signed distances from the input.
     *
     * Note that in some circumstances, this minimum value deviates from
     * what the actual distance function should be, so it may be desirable to do a level set reinitialization after a
     * series of unions.
     *
     * @param   rleFirst                The first rle_level_set operand.
     * @param   rleSecond               The second rle_level_set operand.
     * @param   channelBlendVoxelDistance    The distance, specified in units of voxels, over which the named channel
     * values are blended. A value of 1 or so is usually reasonable.
     */
    void csg_union( const rle_level_set& rleFirst, const rle_level_set& rleSecond,
                    float channelBlendVoxelDistance = 1.5f );

    /**
     * This function sets the current rle level set to the CSG intersection of A and B, two operands.
     *
     * @param   rleFirst                The first rle_level_set operand.
     * @param   rleSecond               The second rle_level_set operand.
     * @param   channelBlendVoxelDistance    The distance, specified in units of voxels, over which the named channel
     * values are blended. A value of 1 or so is usually reasonable.
     */
    void csg_intersect( const rle_level_set& rleFirst, const rle_level_set& rleSecond,
                        float channelBlendVoxelDistance = 1.5f );

    /**
     * This function sets the current rle level set to the CSG subtraction of A and B, two operands.
     * Note that this will not use channels from B in any way. All new voxels in the newly created level set that do not
     * have channel data in A will be set to zero.
     *
     * @param   rleFirst                The first rle_level_set operand.
     * @param   rleSecond               The second rle_level_set operand.
     * @param   channelBlendVoxelDistance    The distance, specified in units of voxels, over which the named channel
     * values are blended. A value of 1 or so is usually reasonable.
     */
    void csg_subtract( const rle_level_set& rleFirst, const rle_level_set& rleSecond,
                       float channelBlendVoxelDistance = 1.5f );

    /**
     * This function linearly interpolates the two rle level sets, as (1-alpha) * rleFirst + alpha * rleSecond.
     *
     * @param   rleFirst                The first rle_level_set operand.
     * @param   rleSecond               The second rle_level_set operand.
     * @param	alpha					The fraction between rleFirst and rleSecond (only values from 0
     * to 1, where 0 returns rleFirst and 1 returns rleSecond).
     *
     */
    void linear_interpolate( const rle_level_set& rleFirst, const rle_level_set& rleSecond, float alpha );

    /**
     * This function linearly interpolates the two rle level sets, as (1-alpha) * rleFirst + alpha * rleSecond, where
     * alpha for each blend is specified by a value in the blendAlphaChannel of the first level set
     *
     * @param   rleFirst                The first rle_level_set operand.
     * @param   rleSecond               The second rle_level_set operand.
     * @param	blendAlphaChannel					The float blend alpha channel from the first
     * level set
     *
     */
    void linear_interpolate( const rle_level_set& rleFirst, const rle_level_set& rleSecond,
                             const frantic::tstring& blendAlphaChannel );

    /**
     * This function resamples a given rle level set into the current voxel coordinate system.
     *
     * @note It is buggy, See Comment Notes!
     *
     * @param	rleSource						The rle_level_set to be resampled.
     * @param	transformNoScaleSourceToCurrent	The rotation and translation from rleSource to the destination rle level
     * set.
     */
    void resample( rle_level_set& rleSource, frantic::graphics::transform4f transformNoScaleSourceToCurrent );

    ////////////////////////////////////////
    // Functions which modify an rle level set in place
    ////////////////////////////////////////

    /**
     * This function permutes the axes of the rle level set.
     *
     * @param  axisPermutation  The permutation for the axes.
     */
    void apply_axis_permutation( const frantic::graphics::vector3& axisPermutation );

    /**
     * This function dilates the defined region of the level set by the specified number of voxels.  It creates a byte
     * channel of the specified name which has value 1 for preexisting defined voxels and 0 for newly defined voxels.
     *
     * @param  dilationVoxels  This is the number of voxels by which to dilate.  It should be 1 or greater.
     * @param  populatedChannelName  The name of the channel to create with 'populated' data.  Set it to "" to not
     * create one.
     */
    void dilate_defined_voxels( int dilationVoxels, const frantic::tstring& populatedChannelName );

    /**
     * This function dilates the defined region of the level set by the specified number of voxels.  It creates a byte
     * channel of the specified name which has value 1 for preexisting defined voxels and 0 for newly defined voxels. It
     * also bounds the result within the provided voxel bounds.
     *
     * @param  dilationVoxels  This is the number of voxels by which to dilate.  It should be 1 or greater.
     * @param  voxelBounds  The defined region will be limited to these bounds.
     * @param  populatedChannelName  The name of the channel to create with 'populated' data.  Set it to "" to not
     * create one.
     */
    void bounded_dilate_defined_voxels( int dilationVoxels, const frantic::graphics::boundbox3& voxelBounds,
                                        const frantic::tstring& populatedChannelName );

    /**
     * This function dilates the defined region of the level set by the values specified in dilationBox.  It creates a
     * byte channel of the specified name which has value 1 for preexisting defined voxels and 0 for newly defined
     * voxels.
     *
     * The dilationBox dilates by boundbox3(-xDilate, +xDilate, -yDilate, +yDilate, -zDilate, +zDilate)
     *
     * @param  dilationBox  This specifies the number of voxels by which to dilate in each of 6 directions.
     * @param  populatedChannelName  The name of the channel to create with 'populated' data.  Set it to "" to not
     * create one.
     */
    void dilate_defined_voxels( const frantic::graphics::boundbox3& dilationBox,
                                const frantic::tstring& populatedChannelName );

  private:
    /**
     * This function is private, and used by dilate_defined_voxels to extend data channels from the old outer bounds.
     *
     * @param  originalOuterBounds  The outer bounds of the rle_index_spec before dilation.
     */
    void extend_channels_from_original_outer_bounds( const frantic::graphics::boundbox3& originalOuterBounds );

  public:
    /**
     * This function trims the level set so that all voxels with a value of 0 in the populated channel become undefined.
     *
     * @param  populatedChannelName  The name of the 'populated' channel to use for trimming.
     */
    void trim_to_populated( const frantic::tstring& populatedChannelName );

    /**
     * This function trims the level set so the outer bounds of the rle_index_spec are intersected with the trimBounds.
     *
     * @param  trimBounds  The bounds for trimming.
     */
    void trim_to_bounds( const frantic::graphics::boundbox3& trimBounds );

    /**
     * This function computes the upwind gradient of a channel. The created gradient channel will have 3 times
     * the arity of the input channel, and if the input channel's arity is not 1, the layout of the gradient will be
     * [all d(data)/dx components], [all d(data)/dy components], [all d(data)/dz components].
     *
     * @param  channelName  The name of the channel for which to compute the gradient.
     * @param  gradientChannelToCreate  The name of the gradient channel to create.
     * @param  populatedChannelToCreate  The populated channel to create.  A voxel is marked as populated (value 1) if
     *                                   all three partial derivatives were estimated.  Pass in an empty string to avoid
     * creating this channel.
     * @param  upwindDir  Which direction to consider "upwind."  +1 means towards higher values, -1 means towards lower
     * values.
     */
    void compute_upwind_gradient( const frantic::tstring& channelName, const frantic::tstring& gradientChannelToCreate,
                                  const frantic::tstring& populatedChannelToCreate, int upwindDir = +1 );

    ////////////////////////////////////////
    // Functions which are used for testing or debugging
    ////////////////////////////////////////

    // This returns the run at the given YZ coordinate
    void getYZ( const rle_level_set& rle, int y, int z );

  private:
    /**
     * This function is private and is used by intersect_ray to find times along the provided ray that are in defined
     * regions. These times are located on the same side of the surface as tMin and tMax, therefore they have the same
     * surface crossing. Since they are defined points, their signed distance values can be used to accurately locate
     * the intersection of the ray and the surface.
     *
     * @param  ray  The ray.
     * @param  tMin  The start time.
     * @param  tMax  The final time.
     * @param  numBisectionSteps  The maximum number of steps to use when doing the binary search for defined regions.
     * @param  tMinOut  A time on the same side of the surface as tMin that is in a defined region.
     * @param  tMaxOut  A time on the same side of the surface as tMax that is in a defined region.
     */
    bool find_defined_times_along_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax,
                                       unsigned int numBisectionSteps, double& tMinOut, double& tMaxOut );

  public:
    /**
     * This function determines the intersection of a ray with the surface of the level set. The return value is a bool
     * of whether there was an intersection.
     *
     * @param  ray  The ray, in world space.
     * @param  tMin  The start time.
     * @param  tMax  The final time.
     * @param  tOut  The time of the surface and ray intersection.
     * @param  outNormal  The surface normal at the intersection point.
     */
    bool intersect_ray( const frantic::graphics::ray3f& ray, double tMin, double tMax, double& tOut,
                        frantic::graphics::vector3f& outNormal );
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic

// Include the header that has the implementations of some template functions
#include <frantic/volumetrics/levelset/rle_level_set_impl.hpp>
