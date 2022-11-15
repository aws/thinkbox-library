// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <map>

#include <frantic/fluids/rle_staggered_vel_accessor.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec_block_iterator_x.hpp>
#include <frantic/volumetrics/levelset/rle_named_channels.hpp>
#include <frantic/volumetrics/voxel_coord_system.hpp>

namespace frantic {
namespace fluids {

/**
 * A rle voxel field is a run-length encoded voxel field that contains any number of data channels. Data channels
 * are generally sampled at voxel centers though it is possible to store data on the voxel faces. In this case
 * each voxel stores 3 samples, while the 3 neighbouring voxels store the remaining face samples
 */
class rle_voxel_field {

    frantic::volumetrics::levelset::rle_index_spec m_rleIndex;
    frantic::volumetrics::voxel_coord_system m_voxelCoordSystem;

    // Generic named channels
    std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel> m_namedChannels;

  public:
    rle_voxel_field() {}

    rle_voxel_field( const frantic::volumetrics::voxel_coord_system& vcs )
        : m_voxelCoordSystem( vcs ) {}

    rle_voxel_field( const frantic::volumetrics::voxel_coord_system& vcs,
                     const frantic::volumetrics::levelset::rle_index_spec& ris )
        : m_rleIndex( ris )
        , m_voxelCoordSystem( vcs )

    {}

    rle_voxel_field( const rle_voxel_field& rhs )
        : m_rleIndex( rhs.m_rleIndex )
        , m_voxelCoordSystem( rhs.m_voxelCoordSystem )
        , m_namedChannels( rhs.m_namedChannels ) {
        // FF_LOG(debug) << "Constructed RLE_Voxel Field Successfully, size " << size() << std::endl;
    }

    /**
     * Swaps the current voxel field with the provided one.
     */
    void swap( rle_voxel_field& rhs ) {
        m_rleIndex.swap( rhs.m_rleIndex );
        m_voxelCoordSystem.swap( rhs.m_voxelCoordSystem );
        m_namedChannels.swap( rhs.m_namedChannels );
    }

    void set_with_swap( frantic::volumetrics::voxel_coord_system& vcs,
                        frantic::volumetrics::levelset::rle_index_spec& ris ) {
        m_rleIndex.swap( ris );
        m_voxelCoordSystem.swap( vcs );
        m_namedChannels.clear();
    }

    /**
     * Swaps an rle_level_set for an rle_voxel_field. Signed distance channel is handled correctly.
     * This is provided mainly as a convienient way to convert to and from level sets and voxel fields.
     * If no "SignedDistance" channel exists in the voxel field, the level set signed distances will be all undefined
     * outside distances.
     * @note  The resulting rle_level_set will have unchanged levelset-specific data such as world-space undefined
     * inside and outside distances, etc.
     *
     * @param levelset  The level set to swap with
     */
    void swap_with_rle_level_set( frantic::volumetrics::levelset::rle_level_set& levelSet );

    void dump( std::ostream& out ) {
        out << "DUMPING RLE VOXEL FIELD\n";
        out << "Defined voxel count: " << m_rleIndex.data_size() << "\n";
        out << "Outer bounds: " << m_rleIndex.outer_bounds() << "\n";
        // m_rleIndex.dump(out);

        out << "Voxel Coordinate System: " << m_voxelCoordSystem << "\n";
        out << "\n";
        out << "Number of named channels: " << m_namedChannels.size() << "\n";
        int index = 0;
        for( std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::const_iterator i =
                 m_namedChannels.begin();
             i != m_namedChannels.end(); ++i ) {
            out << "Channel " << index++ << ", \"" << frantic::strings::to_string( i->first ) << "\"\n";
            out << "Data type: " << frantic::strings::to_string( i->second.type_str() ) << "\n";
            out << "Channel size: " << i->second.size() << "\n";
            // out << "Data: ";
            // channels::channel_data_type_print( out, ",", i->second.arity() * i->second.size(), i->second.data_type(),
            // i->second.data() ); out << "\n";
        }
        out << "FINISHED DUMPING RLE LEVEL SET" << std::endl;
    }

    /**
     * Cleans up the voxel field's data channels and deletes the current rle_index_spec
     */
    void clear() {
        m_rleIndex.clear();
        m_namedChannels.clear();
    }

    /**
     * Set the voxel field to empty, given the exterior region code
     */
    void set_to_empty( int exteriorRegionCode = -1 ) {
        m_rleIndex.clear();
        m_rleIndex.m_exteriorRegionCode = exteriorRegionCode;

        // Clear all the channels to zero-size
        for( std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::iterator i =
                 m_namedChannels.begin();
             i != m_namedChannels.end(); ++i ) {
            i->second.m_data.clear();
        }
    }

    /**
     * Provide access to the voxel coordinate system.
     */
    const frantic::volumetrics::voxel_coord_system& get_voxel_coord_system() const { return m_voxelCoordSystem; }

    /**
     * Provide access to run length encoding index specification.
     */
    const frantic::volumetrics::levelset::rle_index_spec& get_rle_index_spec() const { return m_rleIndex; }

    /**
     * Returns the outer bounds in world space of the field
     */
    frantic::graphics::boundbox3f world_outer_bounds() const {
        return m_voxelCoordSystem.get_world_bounds( m_rleIndex.outer_bounds() );
    }

    /**
     * Returns true if there are not any defined voxels in the voxel field.
     */
    bool empty() const { return m_rleIndex.data_size() == 0; }

    /**
     * Returns the count of defined voxels. Equivalently the size of the data array.
     */
    size_t size() const { return m_rleIndex.data_size(); }

    ////////////////////////////////////////
    // Named channel methods
    ////////////////////////////////////////

    /**
     * Tests the voxel field for a certain channel
     * @param name the channel to look for
     * @returns true if one of the voxel field's named channels has the given name
     */
    bool has_channel( const frantic::tstring& name ) const {
        return m_namedChannels.find( name ) != m_namedChannels.end();
    }

    /**
     * Erases a channel from the voxel field
     * @param name the channel to erase
     */
    void erase_channel( const frantic::tstring& name ) { m_namedChannels.erase( name ); }

    /**
     * Provides a specialized staggered velocity accessor.
     * @param adj - a rle_index_spec adjancency for the current index spec. This is not created behind the scenes so
     * that the caller can reuse the adjacency structure if possible.
     * @throws runtime_error if the voxel field does not have a "StaggeredVelocity" data channel
     */
    const_rle_staggered_vel_accessor get_staggered_vel_accessor() const {
        std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::const_iterator i =
            m_namedChannels.find( _T("StaggeredVelocity") );
        if( i != m_namedChannels.end() ) {
            return const_rle_staggered_vel_accessor( i->second.get_general_accessor(),
                                                     m_rleIndex.get_cached_adjacency() );
        } else {
            throw std::runtime_error( "rle_voxel_field.get_channel_accessor: Tried to retrieve an accessor for channel "
                                      "\"StaggeredVelocity\", but no such channel exists in this rle_level_set." );
        }
    }

    /**
     * Provides a specialized staggered velocity accessor.
     * @param adj - a rle_index_spec adjancency for the current index spec. This is not created behind the scenes so
     * that the caller can reuse the adjacency structure if possible.
     * @throws runtime_error if the voxel field does not have a "StaggeredVelocity" data channel
     */
    rle_staggered_vel_accessor get_staggered_vel_accessor() {

        std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::iterator i =
            m_namedChannels.find( _T("StaggeredVelocity") );
        if( i != m_namedChannels.end() ) {
            return rle_staggered_vel_accessor( i->second.get_general_accessor(), m_rleIndex.get_cached_adjacency() );
        } else {
            throw std::runtime_error( "rle_voxel_field.get_channel_accessor: Tried to retrieve an accessor for channel "
                                      "\"StaggeredVelocity\", but no such channel exists in this rle_level_set." );
        }
    }

    /**
     * This retrieves the names of all the channels.  It adds the names to the existing vector, so if you want
     * just the names this function provides, make sure to initialize the output vector to empty yourself.
     *
     * @param  outNames  A vector into which all the names get added.
     */
    void get_channel_names( std::vector<frantic::tstring>& outNames ) const {
        outNames.reserve( outNames.size() + m_namedChannels.size() );
        for( std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::const_iterator i =
                 m_namedChannels.begin();
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
     * @param  name  The name of the channel to add.
     */
    template <class DataType>
    void add_channel( const frantic::tstring& name ) {
        add_channel( name, frantic::channels::channel_data_type_traits<DataType>::arity(),
                     frantic::channels::channel_data_type_traits<DataType>::data_type() );
    }

    /**
     * This adds a channel to the rle level set, of the given name, data type and arity.
     *
     * It doesn't initialize the data in the named channel, so if you want all zeros or something like
     * that you have to initialize the data yourself.
     *
     * @param  channelName  The name of the channel to add.
     * @param  arity  The arity of the channel.
     * @param  dataType  The data type of the channel.
     */
    void add_channel( const frantic::tstring& channelName, std::size_t arity, frantic::channels::data_type_t dataType );

    /**
     * This sets the specified channel to all zeros.
     *
     * @param  channelName  The name of the channel to set to all zeros.
     */
    void zero_channel( const frantic::tstring& channelName );

    /**
     * This function duplicates the channel named sourceChannelName, copying it into a channel named destChannelName.
     *
     * @param  destChannelName  Where the duplicated channel will end up.
     * @param  sourceChannelName  The source channel to duplicate.
     */
    void duplicate_channel( const frantic::tstring& destChannelName, const frantic::tstring& sourceChannelName );

    /**
     * This function gets the maximum L2-norm of a value in the given channel.
     *
     * @param  channelName  The name of the channel to set to all zeros.
     */
    float get_channel_max_norm( const frantic::tstring& channelName ) const;

    /**
     * This function copies the specified channels from the input voxel field.  Any defined voxels in the destination
     * which are not in the input are set to have channel value 0.  If requested, an uint8 named channel called
     * populatedChannelToCreate is created, which indicates where values were set.  When a 'Populated' channel already
     * exists, and is requested by the caller, the data within the existing channel gets overwritten.
     *
     * @param  inputRVF  The voxel field from which to copy the channels.
     * @param  channelsToCopy  A list of channel names which should be copied.
     * @param  populatedChannelToCreate  If non-empty, the function adds a new uint8 named channel with the given name,
     * as a 'Populated' channel.
     */
    void copy_channels( const rle_voxel_field& inputRVF, const std::vector<frantic::tstring>& channelsToCopy,
                        const frantic::tstring& populatedChannelToCreate = _T("") );

    void copy_channels( const frantic::volumetrics::levelset::rle_level_set& inputRLS,
                        const std::vector<frantic::tstring>& channelsToCopy,
                        const frantic::tstring& populatedChannelToCreate = _T("") );

    frantic::volumetrics::levelset::rle_channel_general_accessor
    get_channel_general_accessor( const frantic::tstring& name ) {
        std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::iterator i =
            m_namedChannels.find( name );
        if( i != m_namedChannels.end() ) {
            return i->second.get_general_accessor();
        } else {
            throw std::runtime_error(
                "rle_voxel_field.get_channel_general_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this rle_voxel_field." );
        }
    }

    frantic::volumetrics::levelset::const_rle_channel_general_accessor
    get_channel_general_accessor( const frantic::tstring& name ) const {
        std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::const_iterator i =
            m_namedChannels.find( name );
        if( i != m_namedChannels.end() ) {
            return i->second.get_general_accessor();
        } else {
            throw std::runtime_error(
                "rle_voxel_field.get_channel_general_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this rle_voxel_field." );
        }
    }

    template <class DataType>
    frantic::volumetrics::levelset::rle_channel_accessor<DataType>
    get_channel_accessor( const frantic::tstring& name ) {
        std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::iterator i =
            m_namedChannels.find( name );
        if( i != m_namedChannels.end() ) {
            return i->second.get_accessor<DataType>();
        } else {
            throw std::runtime_error(
                "rle_voxel_field.get_channel_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this rle_voxel_field." );
        }
    }

    template <class DataType>
    frantic::volumetrics::levelset::const_rle_channel_accessor<DataType>
    get_channel_accessor( const frantic::tstring& name ) const {
        std::map<frantic::tstring, frantic::volumetrics::levelset::rle_channel>::const_iterator i =
            m_namedChannels.find( name );
        if( i != m_namedChannels.end() ) {
            return i->second.get_accessor<DataType>();
        } else {
            throw std::runtime_error(
                "rle_voxel_field.get_channel_accessor: Tried to retrieve an accessor for channel \"" +
                frantic::strings::to_string( name ) + "\", but no such channel exists in this rle_level_set." );
        }
    }

    void linear_interpolate( const rle_voxel_field& rleFirst, const rle_voxel_field& rleSecond, float alpha );

    /**
     * Validates the consistency of this rle_voxel_field.
     *
     * @param  out  The output stream where errors are printed.
     * @return  True if consistent, false otherwise.
     */
    bool check_consistency( std::ostream& out ) const;

    /**
     * Converts a staggered velocity channel to a centered channel in the provided rle level set
     *
     * @param destLevelSet destination level set for the centered velocity channel
     * @param staggeredChannelName  the name of the staggered velocity channel
     * @param centeredChannelName the name of the centered velocity channel
     */
    void convert_staggered_velocity_to_centered( frantic::volumetrics::levelset::rle_level_set& destLevelSet,
                                                 const frantic::tstring& staggeredChannelName,
                                                 const frantic::tstring& centeredChannelName ) const;

    /**
     * Converts a staggered velocity channel to a centered channel in the voxel field
     *
     * @param staggeredChannelName  the name of the staggered velocity channel
     * @param centeredChannelName the name of the centered velocity channel
     */
    void convert_staggered_velocity_to_centered( const frantic::tstring& staggeredChannelName,
                                                 const frantic::tstring& centeredChannelName );

    /**
     * Converts a centered velocity channel in a source rle_level_set to a staggered channel in the voxel field
     *
     * @param sourceLevelSet the level set that contains the centered velocity channel
     * @param centeredChannelName the name of the centered velocity channel
     * @param staggeredChannelName  the name of the staggered velocity channel
     * @param dataIndexMappingChannel the name of the channel in the level set that maps to the corresponding data
     * indices in the voxel field. If not provided a temporary mapping will be created and destroyed
     */
    void convert_centered_velocity_to_staggered( const frantic::volumetrics::levelset::rle_level_set& sourceLevelSet,
                                                 const frantic::tstring& centeredChannelName,
                                                 const frantic::tstring& staggeredChannelName,
                                                 const frantic::tstring& dataIndexMappingChannel );

    void trim_to_bounds( const frantic::graphics::boundbox3& trimVoxelBounds );

    /**
     * This function switches the rle_index_spec being used by this rle_voxel_field.  It sets any new named channel
     * values to 0.
     *
     * This method is copied and modified based on: void rle_level_set::switch_rle_index_spec_with_swap( rle_index_spec&
     * ris, const std::string& populatedChannelToCreate )
     *
     * @param ris  The rle_index_spec to switch to.
     */
    void switch_rle_index_spec_with_swap( frantic::volumetrics::levelset::rle_index_spec& ris,
                                          const frantic::tstring& populatedChannelToCreate = _T("") );

    /**
     * This function permutes the axes of the rle voxel field.
     *
     * @param  axisPermutation  The permutation for the axes.
     */
    void apply_axis_permutation( const frantic::graphics::vector3& axisPermutation );
};

} // namespace fluids
} // namespace frantic
