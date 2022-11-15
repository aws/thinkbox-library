// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/fluids/rle_voxel_field.hpp>
#include <frantic/graphics/boundbox3.hpp>
#include <frantic/graphics/vector3.hpp>

namespace frantic {
namespace fluids {

enum BOUNDARY_CONDITION {
    INVALID = 0,
    NONE = 1,
    DIRICHLET = 2,
    NEUMANN = 4,

};

namespace face_masks {
static const boost::uint32_t X_POS = 15;
static const boost::uint32_t X_NEG = X_POS << 4;
static const boost::uint32_t Y_POS = X_NEG << 4;
static const boost::uint32_t Y_NEG = Y_POS << 4;
static const boost::uint32_t Z_POS = Y_NEG << 4;
static const boost::uint32_t Z_NEG = Z_POS << 4;

inline boost::uint32_t get_shift( boost::uint32_t i ) { return (boost::uint32_t)i * 4; }
inline boost::uint32_t get_mask( boost::uint32_t i ) { return X_POS << i * 4; }

} // namespace face_masks

typedef bool ( *boundary_policy_function_t )( const frantic::graphics::boundbox3& bounds,
                                              const frantic::graphics::vector3& voxelLookup, int face,
                                              boost::uint32_t& outValue );

bool standard_free_boundary_policy( const frantic::graphics::boundbox3& bounds,
                                    const frantic::graphics::vector3& voxelLookup, int face,
                                    boost::uint32_t& outValue );
bool standard_solid_boundary_policy( const frantic::graphics::boundbox3& bounds,
                                     const frantic::graphics::vector3& voxelLookup, int face,
                                     boost::uint32_t& outValue );

void apply_occlusion_boundary_states( const frantic::volumetrics::levelset::rle_index_spec& occlusionRIS,
                                      const float* signedDistanceChannel, boost::int32_t* indexMappingChannel,
                                      boost::uint32_t* faceStateChannel );

void apply_occlusion_boundary( rle_voxel_field& velocityField,
                               frantic::volumetrics::levelset::rle_level_set& occlusionLevelSet,
                               const frantic::tstring& indexMappingChannel = _T("") );

/**
 *  Set the staggered velocity of faces which are inside or which cross
 * the occlusion surface to the occlusion's velocity.
 *
 * @param[in,out] velocityField the velocity field in which to set the
 *      staggered velocities according to the occlusion velocity.
 *      This field must contain a float32[3]
 *      staggered velocity channel, and optionally a uint8 staggered
 *      populated channel.
 * @param staggeredVelocityChannelName the name of the float32[3]
 *      staggered velocity channel in velocityField which shall be
 *      set according to the occlusion velocity.
 * @param staggeredPopulatedChannelName if this is not an empty
 *      string, then this is the name of a uint8 staggered populated
 *      channel which shall be set to 1 for the faces which have
 *      velocities set by this function.
 * @param solidLS the occlusion level set, which may have a
 *      float32[3] "Velocity" channel which specifies the occlusion's
 *      velocity.  Otherwise it is assumed that the occlusion is
 *      stationary.
 * @param velIndexMapChannelName is this is not an empty string, then
 *      this gives the name of a int32 channel in solidLS which maps
 *      the occlusion voxels to the velocity field voxels.  If this
 *      is empty, then such a mapping is created internally by the
 *      function.
 */
void apply_occlusion_boundary_velocity( rle_voxel_field& velocityField,
                                        const frantic::tstring& staggeredVelocityChannelName,
                                        const frantic::tstring& staggeredPopulatedChannelName,
                                        const frantic::volumetrics::levelset::rle_level_set& solidLS,
                                        // const float isosurfaceVoxelDistance,
                                        const frantic::tstring& velIndexMapChannelName = _T("") );

/**
 *  Set the velocity field near the occlusion corresponding to the
 * occlusion boundary condition.
 *
 *  If maintainFluidSeparationChannelName is an empty string, then voxels
 * touching the occlusion are assigned the occlusion's velocity.
 *
 *  Otherwise maintainFluidSeparationChannelName is the name of a uint8
 * channel in the occlusionLS.  Anywhere the separation channel is 1,
 * voxels within outsideVoxelDistance of the occlusion are assigned using
 * the following procedure:
 *   If the velocity is directed away from the occlusion, then the
 *   velocity is kept.
 *   Otherwise, the normal component of the velocity is assigned from
 *   the occlusion's velocity, while the tangential velocity is kept.
 * Anywhere the separation channel is 0, we use the same procedure as if
 * the separationChannelName is empty: voxels touching the occlusion are
 * assigned the occlusion's velocity.
 *
 *  The constrained velocity procedure is taken from
 * frantic::fluids::constrain_staggered_velocity_field().
 *
 * @param[in,out] velocityField the velocity field in which to set the
 *      staggered velocities according to the occlusion velocity.
 *      This field must contain a float32[3] staggered velocity channel.
 * @param staggeredVelocityChannelName the name of the float32[3]
 *      staggered velocity channel in velocityField which shall be
 *      set according to the occlusion velocity.
 * @param occlusionLS the occlusion level set, which may have a
 *      float32[3] "Velocity" channel which specifies the occlusion's
 *      velocity.  Otherwise it is assumed that the occlusion is
 *      stationary.
 * @param outsideVoxelDistance faces with centres closer than this
 *		distance may be altered.  Note that voxels with
 *		separation channel == 0 will always be set if their center
 *		is inside the occlusion.
 * @param maintainFluidSeparationChannelName the name of a uint8 channel
 *		in occlusionLS that specifies whether the velocity should be
 *		constrained to separate from the occlusion.  Voxels set to 0
 *		will have their velocity set to the occlusion's velocity.
 *		Voxels set to 1 will have their velocity set to separate from
 *		the occlusion if it was separaing before, or will have the
 *		normal component of their velocity set to the occlusion's
 *		velocity if it is not separating.
 * @param velIndexMapChannelName is this is not an empty string, then
 *      this gives the name of a int32 channel in solidLS which maps
 *      the occlusion voxels to the velocity field voxels.  If this
 *      is empty, then such a mapping is created internally by the
 *      function.
 */
void apply_constrained_occlusion_boundary_velocity( rle_voxel_field& velocityField,
                                                    const frantic::tstring& staggeredVelocityChannelName,
                                                    const frantic::volumetrics::levelset::rle_level_set& occlusionLS,
                                                    const float outsideVoxelDistance,
                                                    const frantic::tstring& maintainFluidSeparationChannelName,
                                                    const frantic::tstring& velocityIndexMappingChannelName = _T("") );

void apply_sim_boundary( const frantic::graphics::boundbox3& simBounds, rle_voxel_field& velocityField,
                         BOUNDARY_CONDITION wallCondition, const frantic::tstring& faceStateChannel );

/**
 *  Set the staggered velocity on the sim boundary to zero.
 *
 * If staggeredPopulatedChannelName is not an empty string, then the
 * affected faces will be flagged as populated.  Other faces will keep
 * their original flag.
 *
 * @param voxelBounds the simulation boundary.  Voxels will be modified
 *		  if their voxelCoordinate has a component equal to a component
 *		  in voxelBounds.minimum() or in voxelBounds.maximum().
 * @param[in,out] velocityField the field that holds the staggered
 *		  velocity channel to modify.
 * @param staggeredVelocityChannelName the name of a float32[3] staggered
 *        velocity channel in velocityField.
 * @param staggeredPopulatedChannelName the name of an existing uint8
 *		  channel in velocityField.  The bits of this channel are set to
 *		  indicate that the corresponding face in the staggeredVelocity
 *		  channel is populated.  0x01 corresponds to X, 0x02 to Y, and
 *		  0x04 to Z.  This can be an empty string if you do not want
 *		  to modify a populated channel.
 */
void zero_sim_boundary_velocity( const frantic::graphics::boundbox3& voxelBounds,
                                 frantic::fluids::rle_voxel_field& velocityField,
                                 const frantic::tstring& staggeredVelocityChannelName,
                                 const frantic::tstring& staggeredPopulatedChannelName = _T("") );

void apply_fluid_boundary_states( const frantic::volumetrics::levelset::rle_index_spec& fluidRIS,
                                  const float* signedDistanceChannel, boost::int32_t* indexMappingChannel,
                                  boost::uint32_t* faceStateChannel );

void set_free_surface_condition( rle_voxel_field& velocityField,
                                 frantic::volumetrics::levelset::rle_level_set& fluidLevelSet,
                                 const frantic::tstring& faceStateChannel,
                                 const frantic::tstring& indexMappingChannel = _T("") );

/**
 * Sets the face boundary value bits in the uint32. This will reset any previous boundary values.
 *
 *@param value the unsigned int2 that contains face boundary values
 *@param face the face index requested
 *@param boundaryValue the face boundary value
 *@returns the unsigned int2 that contains face boundary values with the bits set
 */
inline boost::uint32_t set_face_boundary( boost::uint32_t value, boost::uint32_t face, boost::uint32_t boundaryValue ) {
    if( face > 5 )
        throw std::runtime_error( "set_face_boundary() - face index " + boost::lexical_cast<std::string>( face ) +
                                  "  provided is not a valid index [0,6)" );

    if( boundaryValue != NONE && boundaryValue != DIRICHLET && boundaryValue != NEUMANN && boundaryValue != INVALID )
        throw std::runtime_error(
            "set_face_boundary() - boundary value " + boost::lexical_cast<std::string>( boundaryValue ) +
            " is not valid, it should be one these values, NEUWMANN=" + boost::lexical_cast<std::string>( NEUMANN ) +
            " DIRICHLET=" + boost::lexical_cast<std::string>( DIRICHLET ) + " NONE=" +
            boost::lexical_cast<std::string>( NONE ) + " INVALID=" + boost::lexical_cast<std::string>( INVALID ) );
    // extract the flag value
    // uint32_t flag = value&face_masks::get_mask(face) >> face_masks::get_shift(face);

    // reset the flag for this face
    value = value & ~face_masks::get_mask( face );

    // TODO we can remove this after some more testing
    // if( (boundaryValue << face_masks::get_shift(face) ) & flag !=0 ) // we have overlapping bits and somethings
    // ismessed 	throw std::runtime_error("bits fucked for boundaryValue= " + lexical_cast<string>(boundaryValue) +
    // " value
    //" + lexical_cast<string>(value));

    return ( boundaryValue << face_masks::get_shift( face ) ) | value;
}

/**
 * Extracts the face boundary value from the containing uint32
 *
 *@param value the unsigned int2 that contains face boundary values
 *@param face the face index requested
 *@returns the face boundary value
 */
inline boost::uint32_t get_face_boundary( boost::uint32_t value, boost::uint32_t face ) {
    if( face > 5 )
        throw std::runtime_error( "face index " + boost::lexical_cast<std::string>( face ) +
                                  "  provided is not a valid index [0,6)" );

    return ( value & face_masks::get_mask( face ) ) >> face_masks::get_shift( face );
}

} // namespace fluids
} // namespace frantic
