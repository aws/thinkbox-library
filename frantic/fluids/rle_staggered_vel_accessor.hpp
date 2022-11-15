// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/graphics/raw_byte_buffer.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/volumetrics/levelset/rle_index_spec.hpp>
#include <frantic/volumetrics/levelset/rle_named_channels.hpp>

namespace frantic {
namespace fluids {

// forward declaration for friending
class rle_voxel_field;
class const_rle_staggered_vel_accessor;

// The non-const accessor adds additional methods for writing to the channel
class rle_staggered_vel_accessor : public volumetrics::levelset::rle_channel_general_accessor {

    // rle_index_spec& m_rleIndexSpec;
    const volumetrics::levelset::ris_adjacency* m_adj;

    /**
     * The usual constructor that gets used by the voxel field to create the accessor.
     * @param src the underlying data accessor to use
     * @param adjacnencies the adjancies of the rle index spec being used
     */
    rle_staggered_vel_accessor( rle_channel_general_accessor src,
                                const volumetrics::levelset::ris_adjacency& adjancencies )
        : rle_channel_general_accessor( src )
        , m_adj( &adjancencies ) {}

    // friend class const_rle_staggered_vel_accessor;
    friend class frantic::volumetrics::levelset::rle_channel;
    friend class rle_voxel_field;
    friend class const_rle_staggered_vel_accessor;

  public:
    /**
     * Copy constructor
     * @param rhs the source velocity accessor to copy
     */
    rle_staggered_vel_accessor( const rle_staggered_vel_accessor& rhs )
        : rle_channel_general_accessor( rhs.m_data, rhs.m_arity, rhs.m_dataType )
        , m_adj( rhs.m_adj ) {}

    /**
     * Assignment operator
     * @param rhs the source velocity accessor
     */
    rle_staggered_vel_accessor& operator=( const rle_staggered_vel_accessor& rhs ) {
        m_data = rhs.m_data;
        m_arity = rhs.m_arity;
        m_dataType = rhs.m_dataType;
        m_primitiveSize = rhs.m_primitiveSize;
        m_weightedSumFunction = rhs.m_weightedSumFunction;

        m_adj = rhs.m_adj;

        return *this;
    }

    /**
     * Default constructor that creates an empty accessor
     */
    rle_staggered_vel_accessor()
        : rle_channel_general_accessor( 0, 0, frantic::channels::data_type_int8 )
        , m_adj( 0 ) {}

    ///////////
    // Access to the data
    ///////////

    /**
     * Array like index access to the underlying data
     * @param i the data index
     */
    const frantic::graphics::vector3f& operator[]( std::size_t i ) const {
        return *reinterpret_cast<const frantic::graphics::vector3f*>(
            m_data->ptr_at( i * sizeof( frantic::graphics::vector3f ) ) );
    }

    /**
     * Array like index access to the underlying data
     * @param i the data index
     */
    frantic::graphics::vector3f& operator[]( std::size_t i ) {
        return *reinterpret_cast<frantic::graphics::vector3f*>(
            m_data->ptr_at( i * sizeof( frantic::graphics::vector3f ) ) );
    }

    /**
     * Adds a new element to the accessed channel
     * @param element the data element to add
     */
    void add_element( const frantic::graphics::vector3f& element ) {
        *reinterpret_cast<frantic::graphics::vector3f*>( rle_channel_general_accessor::add_element() ) = element;
    }

    /**
     * Fills the velocity block provided with the face velocity samples (single floats) indexed by a voxel data index.
     *
     * @param i the data index
     * @param velocities a six element block of floats that will hold the retrieved face velocities
     *        the values are set in this order [ x_i-1/2, x_i+1/2, y_i-1/2, y_i+1/2, z_i-1/2, z_i+1/2 ]
     */
    void get_velocities( std::size_t i, float* velocities ) {

        if( !m_adj ) {
            throw std::runtime_error(
                "rle_staggered_vel_accessor - the rle index spec adjacency is not initialized prior to "
                "being used for a velocity block request" );
        }

        // FF_LOG(debug) << "adj ptr=" << m_adj << std::endl;

        // the adjacency entry
        const frantic::volumetrics::levelset::ris_adj_entry& adj = ( *m_adj )[i];

        // FF_LOG(debug) << "grabbing velocities" << std::endl;
        //  grab the 3 face components indexed at the voxel
        for( int j = 0; j < 3; ++j )
            velocities[j * 2] = *reinterpret_cast<const float*>(
                m_data->ptr_at( ( i * sizeof( frantic::graphics::vector3f ) + j * sizeof( float ) ) ) );

        // FF_LOG(debug) << "grabbed local velocities" << std::endl;
        //  grab the remaining three components from the neighbouring cells
        if( adj.x_pos >= 0 ) {
            // FF_LOG(debug) << "\tadj.xpos " << adj.x_pos << std::endl;
            velocities[1] = *reinterpret_cast<const float*>(
                m_data->ptr_at( ( adj.x_pos * sizeof( frantic::graphics::vector3f ) ) ) );
        }

        // FF_LOG(debug) << "grabbed x_pos vel" << std::endl;
        if( adj.y_pos >= 0 )
            velocities[3] = *reinterpret_cast<const float*>(
                m_data->ptr_at( ( adj.y_pos * sizeof( frantic::graphics::vector3f ) + sizeof( float ) ) ) );
        // FF_LOG(debug) << "grabbed y_pos vel" << std::endl;
        if( adj.z_pos >= 0 )
            velocities[5] = *reinterpret_cast<const float*>(
                m_data->ptr_at( ( adj.z_pos * sizeof( frantic::graphics::vector3f ) + 2 * sizeof( float ) ) ) );
        // FF_LOG(debug) << "grabbed z_pos vel" << std::endl;
    }

    /**
     * Calculates the divergence at the indexed voxel.
     * @param i the data index
     * @param voxelLength the cell size over which to calculate the divergence
     */
    float get_divergence( std::size_t i, float voxelLength ) {
        float vel[6] = { 0, 0, 0, 0, 0, 0 };
        get_velocities( i, vel );
        return ( ( vel[1] - vel[0] ) + ( vel[3] - vel[2] ) + ( vel[5] - vel[4] ) ) / voxelLength;
    }
};

// The non-const accessor adds additional methods for writing to the channel
class const_rle_staggered_vel_accessor : public volumetrics::levelset::const_rle_channel_general_accessor {

    // rle_index_spec& m_rleIndexSpec;
    const volumetrics::levelset::ris_adjacency* m_adj;

    /**
     * The usual constructor that gets used by the voxel field to create the accessor.
     * @param src the underlying data accessor to use
     * @param adjacnencies the adjancies of the rle index spec being used
     */
    const_rle_staggered_vel_accessor( const_rle_channel_general_accessor src,
                                      const volumetrics::levelset::ris_adjacency& adjancencies )
        : const_rle_channel_general_accessor( src )
        , m_adj( &adjancencies ) {}

    // friend class const_rle_staggered_vel_accessor;
    friend class frantic::volumetrics::levelset::rle_channel;
    friend class rle_voxel_field;

  public:
    /**
     * Copy constructor
     * @param rhs the source velocity accessor to copy
     */
    const_rle_staggered_vel_accessor( const rle_staggered_vel_accessor& rhs )
        : const_rle_channel_general_accessor( rhs.m_data, rhs.m_arity, rhs.m_dataType )
        , m_adj( rhs.m_adj ) {}

    /**
     * Assignment operator
     * @param rhs the source velocity accessor
     */
    const_rle_staggered_vel_accessor& operator=( const rle_staggered_vel_accessor& rhs ) {
        m_data = rhs.m_data;
        m_arity = rhs.m_arity;
        m_dataType = rhs.m_dataType;
        m_primitiveSize = rhs.m_primitiveSize;
        m_weightedSumFunction = rhs.m_weightedSumFunction;

        m_adj = rhs.m_adj;

        return *this;
    }

    /**
     * Default constructor that creates an empty accessor
     */
    const_rle_staggered_vel_accessor()
        : const_rle_channel_general_accessor( 0, 0, frantic::channels::data_type_int8 )
        , m_adj( 0 ) {}

    ///////////
    // Access to the data
    ///////////

    /**
     * Array like index access to the underlying data
     * @param i the data index
     */
    const frantic::graphics::vector3f& operator[]( std::size_t i ) const {
        return *reinterpret_cast<const frantic::graphics::vector3f*>(
            m_data->ptr_at( i * sizeof( frantic::graphics::vector3f ) ) );
    }

    /**
     * Fills the velocity block provided with the face velocity samples (single floats) indexed by a voxel data index.
     *
     * @param i the data index
     * @param velocities a six element block of floats that will hold the retrieved face velocities
     *        the values are set in this order [ x_i-1/2, x_i+1/2, y_i-1/2, y_i+1/2, z_i-1/2, z_i+1/2 ]
     */
    void get_velocities( std::size_t i, float* velocities ) const {

        if( !m_adj ) {
            throw std::runtime_error(
                "rle_staggered_vel_accessor - the rle index spec adjacency is not initialized prior to "
                "being used for a velocity block request" );
        }

        // FF_LOG(debug) << "adj ptr=" << m_adj << std::endl;

        // the adjacency entry
        const frantic::volumetrics::levelset::ris_adj_entry& adj = ( *m_adj )[i];

        // FF_LOG(debug) << "grabbing velocities" << std::endl;
        //  grab the 3 face components indexed at the voxel
        for( int j = 0; j < 3; ++j )
            velocities[j * 2] = *reinterpret_cast<const float*>(
                m_data->ptr_at( ( i * sizeof( frantic::graphics::vector3f ) + j * sizeof( float ) ) ) );

        // FF_LOG(debug) << "grabbed local velocities" << std::endl;
        //  grab the remaining three components from the neighbouring cells
        if( adj.x_pos >= 0 ) {
            // FF_LOG(debug) << "\tadj.xpos " << adj.x_pos << std::endl;
            velocities[1] = *reinterpret_cast<const float*>(
                m_data->ptr_at( ( adj.x_pos * sizeof( frantic::graphics::vector3f ) ) ) );
        }

        // FF_LOG(debug) << "grabbed x_pos vel" << std::endl;
        if( adj.y_pos >= 0 )
            velocities[3] = *reinterpret_cast<const float*>(
                m_data->ptr_at( ( adj.y_pos * sizeof( frantic::graphics::vector3f ) + sizeof( float ) ) ) );
        // FF_LOG(debug) << "grabbed y_pos vel" << std::endl;
        if( adj.z_pos >= 0 )
            velocities[5] = *reinterpret_cast<const float*>(
                m_data->ptr_at( ( adj.z_pos * sizeof( frantic::graphics::vector3f ) + 2 * sizeof( float ) ) ) );
        // FF_LOG(debug) << "grabbed z_pos vel" << std::endl;
    }

    /**
     * Calculates the divergence at the indexed voxel.
     * @param i the data index
     * @param voxelLength the cell size over which to calculate the divergence
     */
    float get_divergence( std::size_t i, float voxelLength ) {
        float vel[6] = { 0, 0, 0, 0, 0, 0 };
        get_velocities( i, vel );
        return ( ( vel[1] - vel[0] ) + ( vel[3] - vel[2] ) + ( vel[5] - vel[4] ) ) / voxelLength;
    }
};

} // namespace fluids
} // namespace frantic
