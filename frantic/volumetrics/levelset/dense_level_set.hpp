// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <map>

#include <frantic/graphics2d/boundrect2.hpp>
#include <frantic/volumetrics/levelset/dense_level_set_named_channels.hpp>
#include <frantic/volumetrics/voxel_coord_system.hpp>

namespace frantic {
namespace volumetrics {
namespace levelset {

/**
 * The Dense Level Set uses an single block of memory to store values for all the voxels in a box region.
 *
 * @todo  Will have to sort out narrow-band details, for now those considerations are being fully ignored.
 * @todo  Will have to add an exterior region code to the dense level set as well.
 */
class dense_level_set {
    voxel_coord_system m_voxelCoordSystem;

    // These two parameters specify the bounding box of the data array in voxel space
    frantic::graphics::vector3 m_voxelBoundsMin;
    frantic::graphics::size3 m_voxelBoundsSize;

    // The distance channel, treated specially
    std::vector<float> m_distanceData;

    // Generic named channels
    std::map<frantic::tstring, dense_level_set_channel> m_namedChannels;

  public:
    dense_level_set()
        : m_voxelBoundsMin( 0 )
        , m_voxelBoundsSize( 0 ) {}

    dense_level_set( const voxel_coord_system& vcs )
        : m_voxelCoordSystem( vcs )
        , m_voxelBoundsMin( 0 )
        , m_voxelBoundsSize( 0 ) {}

    // Returns true if there are any defined voxels in the level set.
    bool empty() const { return m_distanceData.empty(); }

    // Provide array-based access to the data
    std::size_t size() const { return m_distanceData.size(); }

    frantic::graphics::boundbox3 outer_bounds() const {
        return frantic::graphics::boundbox3( m_voxelBoundsMin, m_voxelBoundsSize );
    }

    frantic::graphics::boundbox3f world_outer_bounds() const {
        return m_voxelCoordSystem.get_world_bounds( outer_bounds() );
    }

    void swap( dense_level_set& rhs ) {
        m_voxelCoordSystem.swap( rhs.m_voxelCoordSystem );
        std::swap( m_voxelBoundsMin, rhs.m_voxelBoundsMin );
        std::swap( m_voxelBoundsSize, rhs.m_voxelBoundsSize );
        m_distanceData.swap( rhs.m_distanceData );
        m_namedChannels.swap( rhs.m_namedChannels );
    }

    void clear() {
        m_voxelBoundsMin.set( 0 );
        m_voxelBoundsSize.set( 0 );
        m_distanceData.clear();
        m_namedChannels.clear();
    }

    // Provide array-based access to the data
    float& operator[]( int index ) { return m_distanceData[index]; }

    // Provide array-based access to the data
    const float& operator[]( int index ) const { return m_distanceData[index]; }

    // Provide an [] operator using voxel coordinate.  Throws an exception when trying to access invalid data.
    float& operator[]( const frantic::graphics::vector3& voxel ) {
        int dataIndex = m_voxelBoundsSize.get_index_nothrow( voxel - m_voxelBoundsMin );
        if( dataIndex >= 0 )
            return m_distanceData[dataIndex];
        else
            throw std::runtime_error( "dense_level_set.operator[]: Tried to access out of bounds voxel coordinate " +
                                      voxel.str() + ", in a level set whose bounding box reaches from " +
                                      m_voxelBoundsMin.str() + " to " +
                                      ( m_voxelBoundsMin + m_voxelBoundsSize ).str() );
    }

    // Provide an [] operator using voxel coordinate.  Throws an exception when trying to access invalid data.
    const float& operator[]( const frantic::graphics::vector3& voxel ) const {
        int dataIndex = m_voxelBoundsSize.get_index_nothrow( voxel - m_voxelBoundsMin );
        if( dataIndex >= 0 )
            return m_distanceData[dataIndex];
        else
            throw std::runtime_error( "dense_level_set.operator[]: Tried to access out of bounds voxel coordinate " +
                                      voxel.str() + ", in a level set whose bounding box reaches from " +
                                      m_voxelBoundsMin.str() + " to " +
                                      ( m_voxelBoundsMin + m_voxelBoundsSize ).str() );
    }

    // Provide an == operator to compare two dense level sets
    bool operator==( const dense_level_set& dls ) {
        if( !( m_voxelCoordSystem == dls.m_voxelCoordSystem ) )
            return false;

        if( m_voxelBoundsMin != dls.m_voxelBoundsMin || m_voxelBoundsSize != dls.m_voxelBoundsSize )
            return false;

        if( !( m_distanceData == dls.m_distanceData ) )
            return false;

        //		if( !( m_namedChannels == dls.m_namedChannels ) )
        //			return false;

        return true;
    }

    voxel_coord_system& get_voxel_coord_system() { return m_voxelCoordSystem; }

    const voxel_coord_system& get_voxel_coord_system() const { return m_voxelCoordSystem; }

    const std::vector<float> get_distance_data() const { return m_distanceData; }

    ////////////////////////////////////////
    // Named channel methods
    ////////////////////////////////////////

    bool has_channel( const frantic::tstring& name ) const {
        return m_namedChannels.find( name ) != m_namedChannels.end();
    }

    void erase_channel( const frantic::tstring& name ) { m_namedChannels.erase( name ); }

    // This retrieves the names of all the channels
    void get_channel_names( std::vector<frantic::tstring>& outNames ) const {
        outNames.reserve( outNames.size() + m_namedChannels.size() );
        for( std::map<frantic::tstring, dense_level_set_channel>::const_iterator i = m_namedChannels.begin();
             i != m_namedChannels.end(); ++i ) {
            outNames.push_back( i->first );
        }
    }

    template <class DataType>
    void add_channel( const frantic::tstring& name ) {
        add_channel( name, frantic::channels::channel_data_type_traits<DataType>::arity(),
                     frantic::channels::channel_data_type_traits<DataType>::data_type() );
    }

    void add_channel( const frantic::tstring& name, std::size_t arity, data_type_t dataType ) {
        std::map<frantic::tstring, dense_level_set_channel>::iterator i = m_namedChannels.find( name );
        if( i == m_namedChannels.end() ) {
            m_namedChannels.insert( std::make_pair( name, dense_level_set_channel( name, arity, dataType ) ) );
            i = m_namedChannels.find( name );
            // Make sure the vertex array count matches that of the rle_index_spec's data size
            i->second.m_data.resize( m_distanceData.size() * i->second.primitive_size() );
        } else {
            // TODO: Should this rather erase the existing channel and replace it with a new one?
            if( i->second.arity() != arity || i->second.data_type() != dataType )
                throw std::runtime_error(
                    "dense_level_set.add_channel: Tried to add channel \"" + frantic::strings::to_string( name ) +
                    "\", with data type " + frantic::strings::to_string( channel_data_type_str( arity, dataType ) ) +
                    ".  This could not be done, because the channel already exists with data type " +
                    frantic::strings::to_string( channel_data_type_str( i->second.arity(), i->second.data_type() ) ) +
                    "." );
        }
    }

    dense_level_set_channel_general_accessor get_channel_general_accessor( const frantic::tstring& name ) {
        std::map<frantic::tstring, dense_level_set_channel>::iterator i = m_namedChannels.find( name );
        if( i != m_namedChannels.end() ) {
            return i->second.get_general_accessor();
        } else {
            throw std::runtime_error(
                "dense_level_set.get_channel_general_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this dense_level_set." );
        }
    }

    const_dense_level_set_channel_general_accessor get_channel_general_accessor( const frantic::tstring& name ) const {
        std::map<frantic::tstring, dense_level_set_channel>::const_iterator i = m_namedChannels.find( name );
        if( i != m_namedChannels.end() ) {
            return i->second.get_general_accessor();
        } else {
            throw std::runtime_error(
                "dense_level_set.get_channel_general_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this dense_level_set." );
        }
    }

    template <class DataType>
    dense_level_set_channel_accessor<DataType> get_channel_accessor( const frantic::tstring& name ) {
        std::map<frantic::tstring, dense_level_set_channel>::iterator i = m_namedChannels.find( name );
        if( i != m_namedChannels.end() ) {
            return i->second.get_accessor<DataType>();
        } else {
            throw std::runtime_error(
                "dense_level_set.get_channel_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this dense_level_set." );
        }
    }

    template <class DataType>
    const_dense_level_set_channel_accessor<DataType> get_channel_accessor( const frantic::tstring& name ) const {
        std::map<frantic::tstring, dense_level_set_channel>::const_iterator i = m_namedChannels.find( name );
        if( i != m_namedChannels.end() ) {
            return i->second.get_accessor<DataType>();
        } else {
            throw std::runtime_error(
                "dense_level_set.get_channel_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this dense_level_set." );
        }
    }

    ////////////////////////////////////////
    // Queries to get data
    ////////////////////////////////////////

    // This function fills a dense xy plane as specified by the bounding rectangle with sample values from the level
    // set.
    void fill_plane( const frantic::graphics2d::boundrect2& voxelXYExtents, int voxelZ,
                     std::vector<float>& outVoxelCornerValues ) const;

    // This function fills a dense integer box with sample values from the level set.
    void fill_box( const frantic::graphics::boundbox3& voxelExtents, std::vector<float>& outVoxelCornerValues ) const;

    /*
     * Dumps the data of the dense_level_set to the provided output stream.
     *
     * @param   out  The ostream the structure is dumped to.
     */
    void dump( std::ostream& out ) const;

    ////////////////////////////////////////
    // Functions which construct an dense level set
    ////////////////////////////////////////

    // This sets the dense level set to the given data, using swap functions to set the data efficiently.  Note that
    // the variables you pass in will have arbitrary values in them after this function is called.
    void set_with_swap( voxel_coord_system& vcs, const graphics::boundbox3& voxelBounds,
                        std::vector<float>& distanceData );

    /**
     * This function overwrites the current dense level set with the CSG union of the two operands.  This is done by
     * taking the minimum value of the signed distances from the input.
     *
     * Note that in some circumstances, this minimum value deviates from
     * what the actual distance function should be, so it may be desirable to do a level set reinitialization after a
     * series of unions.
     *
     * @param   dlsFirst                The first dense_level_set operand.
     * @param   dlsSecond               The second dense_level_set operand.
     * @param   channelBlendDistance    The distance, specified in units of voxels, over which the named channel values
     * are blended. A value of 1 or so is usually reasonable.
     */
    //	void csg_union(const dense_level_set& dlsFirst, const dense_level_set& dlsSecond, float channelBlendDistance
    //= 1.f);

    // This function sets the current dense level set to the common areas of A and B, two previous level sets
    // NOTE: This complements the operands, does a union, then complements the operands again.  This is why they are not
    //       const, though after it's done the operands will be back to the value they started.  This affects any usage
    //       of this function for multithreading purposes.
    //	void csg_intersect(dense_level_set& dlsFirst, dense_level_set& dlsSecond, float channelBlendDistance = 1.f);

    ////////////////////////////////////////
    // Functions which modify an dense level set in place
    ////////////////////////////////////////

    // This function converts an dense level set to its complement
    //	void csg_complement();
};

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
