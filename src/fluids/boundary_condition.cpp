// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/foreach.hpp>

#include <frantic/fluids/boundary_condition.hpp>

#include <frantic/geometry/trimesh3_file_io.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/volumetrics/levelset/rle_defined_and_adj_iterator.hpp>
#include <frantic/volumetrics/levelset/rle_level_set.hpp>
#include <frantic/volumetrics/levelset/rle_trilerp.hpp>

#include <frantic/files/filename_sequence.hpp>

using namespace std;
using namespace boost;

using namespace frantic::graphics;
using namespace frantic::fluids;
using namespace frantic::volumetrics::levelset;

void frantic::fluids::apply_fluid_boundary_states( const frantic::volumetrics::levelset::rle_index_spec& fluidRIS,
                                                   const float* signedDistanceChannel,
                                                   boost::int32_t* indexMappingChannel,
                                                   boost::uint32_t* faceStateChannel ) {
    const ris_adjacency& fluidAdj = fluidRIS.get_cached_adjacency();

    for( size_t i = 0, ie = fluidAdj.data_size(); i != ie; ++i ) {
        boost::int32_t velocityFieldIndex = indexMappingChannel[i];
        if( velocityFieldIndex >= 0 ) {
            float fluidDistance = signedDistanceChannel[i];
            const ris_adj_entry& adj = fluidAdj[i];

            for( int face = 0; face < 6; ++face ) {
                int neighborIndex = adj[face];

                // if the neighbor is outside set the face to the free condition
                if( neighborIndex == -1 ) {
                    faceStateChannel[velocityFieldIndex] =
                        set_face_boundary( faceStateChannel[velocityFieldIndex], face, DIRICHLET );
                } else if( neighborIndex >= 0 ) {
                    const float neighborDistance = signedDistanceChannel[neighborIndex];
                    // int neighborVelocityIndex = indexMappingChannel[neighborIndex];
                    //  if the voxel is outside then all the faces (regardless of in or out) is set to the free
                    //  condition
                    if( fluidDistance > 0.f ) {
                        faceStateChannel[velocityFieldIndex] =
                            set_face_boundary( faceStateChannel[velocityFieldIndex], face, DIRICHLET );
                    }
                    // the voxel is a zero crossing from inside the fluid to the air
                    else if( neighborDistance > 0.f ) {
                        faceStateChannel[velocityFieldIndex] =
                            set_face_boundary( faceStateChannel[velocityFieldIndex], face, DIRICHLET );
                    }
                }
            }
        }
    }
}

void frantic::fluids::set_free_surface_condition( rle_voxel_field& velocityField,
                                                  frantic::volumetrics::levelset::rle_level_set& fluidLevelSet,
                                                  const frantic::tstring& faceStateChannel,
                                                  const frantic::tstring& indexMappingChannel ) {
    if( fluidLevelSet.empty() || velocityField.empty() )
        return;

    const rle_index_spec& fluidRIS = fluidLevelSet.get_rle_index_spec();
    const rle_index_spec& velocityRIS = velocityField.get_rle_index_spec();

    std::vector<boost::int32_t> indexMapping;

    boost::int32_t* indexMappingPtr = 0;

    if( indexMappingChannel.empty() ) {
        indexMapping.resize( fluidRIS.data_size() );
        indexMappingPtr = &indexMapping[0];

        fluidRIS.fill_data_index_map( velocityRIS, indexMappingPtr );
    } else if( !fluidLevelSet.has_channel( indexMappingChannel ) ) {
        throw std::runtime_error( "set_free_surface_condition() - Data Index Mapping Channel \"" +
                                  frantic::strings::to_string( indexMappingChannel ) +
                                  "\" does not exist in the occlusion level set provided" );
    } else {
        rle_channel_accessor<boost::int32_t> mappingAcc =
            fluidLevelSet.get_channel_accessor<boost::int32_t>( indexMappingChannel );
        indexMappingPtr = &mappingAcc[0];
    }

    if( !velocityField.has_channel( faceStateChannel ) ) {
        throw std::runtime_error( "set_free_surface_condition() - voxel field must have a \"" +
                                  frantic::strings::to_string( faceStateChannel ) + "\" channel for this operation" );
    }

    rle_channel_accessor<boost::uint32_t> faceStateAcc =
        velocityField.get_channel_accessor<boost::uint32_t>( faceStateChannel );

    boost::uint32_t* faceStatePtr = reinterpret_cast<boost::uint32_t*>( faceStateAcc.data( 0 ) );

    apply_fluid_boundary_states( fluidRIS, &fluidLevelSet[0], indexMappingPtr, faceStatePtr );
}

void frantic::fluids::apply_occlusion_boundary_states( const rle_index_spec& occlusionRIS,
                                                       const float* signedDistanceChannel,
                                                       boost::int32_t* indexMappingChannel,
                                                       boost::uint32_t* faceStateChannel ) {

    const ris_adjacency& occAdj = occlusionRIS.get_cached_adjacency();

    for( size_t i = 0, ie = occAdj.data_size(); i != ie; ++i ) {

        int velocityFieldIndex = indexMappingChannel[i];
        if( velocityFieldIndex < 0 )
            continue;

        float occlusionDistance = signedDistanceChannel[i];

        const ris_adj_entry& adj = occAdj[i];

        for( int face = 0; face < 6; ++face ) {

            int neighborIndex = adj[face];

            // if the neighbor is inside set the face to invalid
            if( neighborIndex <= -2 ) {
                faceStateChannel[velocityFieldIndex] =
                    set_face_boundary( faceStateChannel[velocityFieldIndex], face, NEUMANN );
            } else if( neighborIndex >= 0 ) {
                const float occlusionNeighborDistance = signedDistanceChannel[neighborIndex];
                // all faces which are inside or cross the occlusion get the NEUMANN condition
                if( occlusionDistance <= 0 || occlusionNeighborDistance <= 0 )
                    faceStateChannel[velocityFieldIndex] =
                        set_face_boundary( faceStateChannel[velocityFieldIndex], face, NEUMANN );
            }
        }
    }
}

void frantic::fluids::apply_occlusion_boundary( rle_voxel_field& velocityField, rle_level_set& occlusionLevelSet,
                                                const frantic::tstring& indexMappingChannel ) {
    const rle_index_spec& occRIS = occlusionLevelSet.get_rle_index_spec();
    const rle_index_spec& velocityRIS = velocityField.get_rle_index_spec();

    if( occlusionLevelSet.empty() )
        return;

    std::vector<boost::int32_t> indexMapping;

    boost::int32_t* indexMappingPtr = 0;

    if( indexMappingChannel.empty() ) {
        indexMapping.resize( occRIS.data_size() );
        indexMappingPtr = &indexMapping[0];

        occRIS.fill_data_index_map( velocityRIS, indexMappingPtr );
    } else if( !occlusionLevelSet.has_channel( indexMappingChannel ) ) {
        throw std::runtime_error( "apply_occlusion_boundary() - Data Index Mapping Channel \"" +
                                  frantic::strings::to_string( indexMappingChannel ) +
                                  "\" does not exist in the occlusion level set provided" );
    } else {
        rle_channel_accessor<boost::int32_t> mappingAcc =
            occlusionLevelSet.get_channel_accessor<boost::int32_t>( indexMappingChannel );
        indexMappingPtr = &mappingAcc[0];
    }

    if( !velocityField.has_channel( _T("FaceState") ) ) {
        throw std::runtime_error(
            "apply_occlusion_boundary() - voxel field must have a \"FaceState\" channel for this operation" );
    }

    if( !velocityField.has_channel( _T("StaggeredVelocity") ) ) {
        throw std::runtime_error(
            "apply_occlusion_boundary() - voxel field must have a \"StaggeredVelocity\" channel for this operation" );
    }

    rle_channel_accessor<boost::uint32_t> faceStateAcc =
        velocityField.get_channel_accessor<boost::uint32_t>( _T("FaceState") );
    rle_channel_accessor<vector3f> staggeredAcc =
        velocityField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    FF_LOG( error ) << "**** (pre occ constraint) staggered velocity channel size=" << staggeredAcc.size() << " ****"
                    << endl;

    boost::uint32_t* faceStatePtr = reinterpret_cast<boost::uint32_t*>( faceStateAcc.data( 0 ) );
    //	vector3f* velocityPtr = reinterpret_cast<vector3f*>(staggeredAcc.data(0));

    FF_LOG( error ) << "**** (pre occ constraint [post ptr referencing ) staggered velocity channel size="
                    << staggeredAcc.size() << " ****" << endl;
    apply_occlusion_boundary_states( occRIS, &occlusionLevelSet[0], indexMappingPtr, faceStatePtr );
    FF_LOG( error ) << "**** (post occ constraint) staggered velocity channel size=" << staggeredAcc.size() << " ****"
                    << endl;
}

void frantic::fluids::apply_occlusion_boundary_velocity( rle_voxel_field& velocityField,
                                                         const frantic::tstring& staggeredVelocityChannelName,
                                                         const frantic::tstring& staggeredPopulatedChannelName,
                                                         const frantic::volumetrics::levelset::rle_level_set& solidLS,
                                                         const frantic::tstring& velIndexMapChannelName ) {
    const rle_index_spec& solidRIS = solidLS.get_rle_index_spec();
    const rle_index_spec& velocityRIS = velocityField.get_rle_index_spec();

    if( !velocityField.has_channel( staggeredVelocityChannelName ) )
        throw std::runtime_error(
            "apply_occlusion_boundary_velocity Error: the velocity field is missing the specified \'" +
            frantic::strings::to_string( staggeredVelocityChannelName ) + "\' staggered velocity channel." );

    rle_channel_accessor<vector3f> velAcc =
        velocityField.get_channel_accessor<vector3f>( staggeredVelocityChannelName );

    const boost::int32_t* velIndexMap = 0;
    std::vector<boost::int32_t> velIndexMapVector;
    if( velIndexMapChannelName.empty() ) {
        velIndexMapVector.resize( solidRIS.data_size() );
        solidRIS.fill_data_index_map( velocityRIS, &velIndexMapVector[0] );
        velIndexMap = &velIndexMapVector[0];
    } else if( !solidLS.has_channel( velIndexMapChannelName ) ) {
        throw std::runtime_error(
            "apply_occlusion_boundary_velocity Error: the occlusion level set is missing the specified \'" +
            frantic::strings::to_string( velIndexMapChannelName ) +
            "\' occlusion to velocity field index map channel." );
    } else {
        const_rle_channel_accessor<boost::int32_t> velIndexMapAcc =
            solidLS.get_channel_accessor<boost::int32_t>( velIndexMapChannelName );
        velIndexMap = &velIndexMapAcc[0];
    }

    const bool solidHasVelocity = solidLS.has_channel( _T("Velocity") );
    const_rle_channel_accessor<vector3f> solidVelAcc;

    if( solidHasVelocity ) {
        solidVelAcc = solidLS.get_channel_accessor<vector3f>( _T("Velocity") );
    }

    const bool hasPopulatedChannel = !staggeredPopulatedChannelName.empty();
    rle_channel_accessor<unsigned char> populatedAcc;

    if( hasPopulatedChannel ) {
        if( !velocityField.has_channel( staggeredPopulatedChannelName ) )
            velocityField.add_channel<unsigned char>( staggeredPopulatedChannelName );

        populatedAcc = velocityField.get_channel_accessor<unsigned char>( staggeredPopulatedChannelName );
    }

    const ris_adjacency& solidAdj = solidRIS.get_cached_adjacency();

    for( rle_defined_iterator i( solidRIS ), ie( solidRIS, true ); i != ie; ++i ) {
        const boost::int32_t solidIndex = i.get_data_index();
        const boost::int32_t velIndex = velIndexMap[solidIndex];
        if( velIndex >= 0 ) {
            const float centerPhi = solidLS[solidIndex];
            const vector3f centerVel = solidHasVelocity ? solidVelAcc[solidIndex] : vector3f( 0 );

            for( int axis = 0; axis < 3; ++axis ) {
                const int face = get_negative_rae_neighbor_index( axis );
                const boost::int32_t adjSolidIndex = solidAdj[solidIndex][face];
                if( adjSolidIndex <= -2 ) {
                    // adj is deep within the occlusion
                    velAcc[velIndex][axis] = centerVel[axis];
                    if( hasPopulatedChannel )
                        populatedAcc[velIndex] |= 0x01 << axis;
                } else if( adjSolidIndex >= 0 ) {
                    const float adjPhi = solidLS[adjSolidIndex];
                    // adj is defined
                    if( centerPhi <= 0 || adjPhi <= 0 ) {
                        const float adjVel = solidHasVelocity ? solidVelAcc[adjSolidIndex][axis] : 0;
                        const float faceVel = 0.5f * ( centerVel[axis] + adjVel );
                        velAcc[velIndex][axis] = faceVel;
                        if( hasPopulatedChannel )
                            populatedAcc[velIndex] |= 0x01 << axis;
                    }
                }
            }
        }
    }
}

void frantic::fluids::apply_constrained_occlusion_boundary_velocity(
    rle_voxel_field& velocityField, const frantic::tstring& staggeredVelocityChannelName,
    const rle_level_set& occlusionLS, const float outsideVoxelDistance,
    const frantic::tstring& maintainFluidSeparationChannelName,
    const frantic::tstring& velocityIndexMappingChannelName ) {
    if( occlusionLS.size() == 0 )
        return;

    typedef boost::uint8_t separation_t;

    const bool defaultFluidSeparation = false;
    const float voxelLength = occlusionLS.get_voxel_coord_system().voxel_length();
    const float isosurfaceDistance = outsideVoxelDistance * voxelLength;
    const float assignOcclusionVelocityIsosurfaceDistance = -0.5f * voxelLength;

    if( staggeredVelocityChannelName.empty() ) {
        throw std::runtime_error(
            "apply_occlusion_boundary_velocity Error: missing required parameter staggered velocity channel name" );
    }
    if( !velocityField.has_channel( staggeredVelocityChannelName ) ) {
        throw std::runtime_error(
            "apply_occlusion_boundary_velocity Error: the velocity field is missing the specified \'" +
            frantic::strings::to_string( staggeredVelocityChannelName ) + "\' staggered velocity channel." );
    } else {
        typedef frantic::channels::channel_data_type_traits<vector3f> required_traits;
        rle_channel_general_accessor tempAcc =
            velocityField.get_channel_general_accessor( staggeredVelocityChannelName );
        if( tempAcc.arity() != required_traits::arity() || tempAcc.data_type() != required_traits::data_type() )
            throw std::runtime_error(
                "apply_occlusion_boundary_velocity Error: the specified \'" +
                frantic::strings::to_string( staggeredVelocityChannelName ) + "\' channel must be of type " +
                frantic::strings::to_string( required_traits::type_str() ) + " but instead it is of type " +
                frantic::strings::to_string( tempAcc.type_str() ) + "." );
    }
    rle_channel_accessor<vector3f> inputStaggeredVelocityAcc =
        velocityField.get_channel_accessor<vector3f>( staggeredVelocityChannelName );
    std::vector<vector3f> outputStaggeredVelocity( inputStaggeredVelocityAcc.size() );
    memcpy( &outputStaggeredVelocity[0], &inputStaggeredVelocityAcc[0],
            inputStaggeredVelocityAcc.size() * sizeof( vector3f ) );

    const separation_t* separationChannel = 0;
    if( !maintainFluidSeparationChannelName.empty() ) {
        typedef frantic::channels::channel_data_type_traits<separation_t> required_traits;
        const_rle_channel_general_accessor tempAcc =
            occlusionLS.get_channel_general_accessor( maintainFluidSeparationChannelName );
        if( tempAcc.arity() != required_traits::arity() || tempAcc.data_type() != required_traits::data_type() )
            throw std::runtime_error( "apply_occlusion_boundary_velocity Error: the specified \'" +
                                      frantic::strings::to_string( maintainFluidSeparationChannelName ) +
                                      "\' channel in the occlusion level set must be of type " +
                                      frantic::strings::to_string( required_traits::type_str() ) +
                                      " but instead it is of type " +
                                      frantic::strings::to_string( tempAcc.type_str() ) + "." );
        else {
            const_rle_channel_accessor<separation_t> separationAcc =
                occlusionLS.get_channel_accessor<separation_t>( maintainFluidSeparationChannelName );
            separationChannel = &separationAcc[0];
        }
    }

    const rle_index_spec& velocityFieldRIS = velocityField.get_rle_index_spec();
    const rle_index_spec& occlusionRIS = occlusionLS.get_rle_index_spec();
    const ris_adjacency& occlusionAdj = occlusionRIS.get_cached_adjacency();

    // temporary storage in case named channel is not allocated
    std::vector<boost::int32_t> velocityIndexMappingVector;
    const boost::int32_t* velocityIndexMap = 0;

    if( velocityIndexMappingChannelName.empty() ) {
        velocityIndexMappingVector.resize( occlusionLS.size() );
        occlusionRIS.fill_data_index_map( velocityFieldRIS, &velocityIndexMappingVector[0] );
        velocityIndexMap = &velocityIndexMappingVector[0];
    } else if( !occlusionLS.has_channel( velocityIndexMappingChannelName ) ) {
        throw std::runtime_error(
            "apply_occlusion_boundary_velocity Error: the occlusion level set is missing the specified \'" +
            frantic::strings::to_string( velocityIndexMappingChannelName ) +
            "\' occlusion to velocity index mapping channel." );
    } else {
        typedef frantic::channels::channel_data_type_traits<boost::int32_t> required_traits;
        const_rle_channel_general_accessor tempAcc =
            occlusionLS.get_channel_general_accessor( velocityIndexMappingChannelName );
        if( tempAcc.arity() != required_traits::arity() || tempAcc.data_type() != required_traits::data_type() )
            throw std::runtime_error( "apply_occlusion_boundary_velocity Error: the specified \'" +
                                      frantic::strings::to_string( velocityIndexMappingChannelName ) +
                                      "\' channel must in the occlusion level set be of type " +
                                      frantic::strings::to_string( required_traits::type_str() ) +
                                      " but instead it is of type " +
                                      frantic::strings::to_string( tempAcc.type_str() ) + "." );
        else {
            const_rle_channel_accessor<boost::int32_t> velocityIndexMapAcc =
                occlusionLS.get_channel_accessor<boost::int32_t>( velocityIndexMappingChannelName );
            velocityIndexMap = &velocityIndexMapAcc[0];
        }
    }

    const frantic::tstring velocityChannelName = _T("Velocity");
    const bool occlusionHasVelocity = occlusionLS.has_channel( velocityChannelName );
    const_rle_channel_accessor<vector3f> occlusionVelocityAcc;
    if( occlusionHasVelocity ) {
        typedef frantic::channels::channel_data_type_traits<vector3f> required_traits;
        const_rle_channel_general_accessor tempAcc = occlusionLS.get_channel_general_accessor( velocityChannelName );
        if( tempAcc.arity() != required_traits::arity() || tempAcc.data_type() != required_traits::data_type() )
            throw std::runtime_error(
                "apply_occlusion_boundary_velocity Error: the \'" + frantic::strings::to_string( velocityChannelName ) +
                "\' channel must in the occlusion level set be of type " +
                frantic::strings::to_string( required_traits::type_str() ) + " but instead it is of type " +
                frantic::strings::to_string( tempAcc.type_str() ) + "." );
        else {
            occlusionVelocityAcc = occlusionLS.get_channel_accessor<vector3f>( velocityChannelName );
        }
    }

    for( rle_defined_iterator i = occlusionRIS.begin(); i != occlusionRIS.end(); ++i ) {
        const boost::int32_t occlusionIndex = i.get_data_index();
        const boost::int32_t velocityIndex = velocityIndexMap[occlusionIndex];

        if( velocityIndex < 0 )
            continue;

        const vector3 voxelCoord = i.get_coord();

        const vector3f initialVelocity = inputStaggeredVelocityAcc[velocityIndex];
        vector3f resultVelocity = initialVelocity;

        for( int axis = 0; axis < 3; ++axis ) {
            const int face = get_negative_rae_neighbor_index( axis );
            const int adjOcclusionIndex = occlusionAdj[occlusionIndex][face];

            if( adjOcclusionIndex <= -2 ) {
                // deep inside occlusion
                resultVelocity[axis] = 0;
            } else if( adjOcclusionIndex >= 0 ) {
                const bool ensureFluidSeparation = separationChannel == 0 ? defaultFluidSeparation
                                                                          : ( separationChannel[occlusionIndex] ||
                                                                              separationChannel[adjOcclusionIndex] );
                if( ensureFluidSeparation ) {
                    const float occlusionDistance =
                        0.5f * ( occlusionLS[occlusionIndex] + occlusionLS[adjOcclusionIndex] );

                    if( occlusionDistance <= isosurfaceDistance ) {
                        vector3f offset;
                        switch( axis ) {
                        case 0:
                            offset = vector3f( 0.f, 0.5f, 0.5f );
                            break;
                        case 1:
                            offset = vector3f( 0.5f, 0.f, 0.5f );
                            break;
                        case 2:
                            offset = vector3f( 0.5f, 0.5f, 0.f );
                            break;
                        }
                        const vector3f faceCoord = voxelCoord + offset;

                        vector3f occlusionVelocityAtFace( 0 );
                        if( occlusionHasVelocity ) {
                            occlusionVelocityAtFace = 0.5f * ( occlusionVelocityAcc[occlusionIndex] +
                                                               occlusionVelocityAcc[adjOcclusionIndex] );
                        }

                        if( occlusionDistance < assignOcclusionVelocityIsosurfaceDistance ) {
                            resultVelocity[axis] = occlusionVelocityAtFace[axis];
                        } else {
                            vector3f occlusionGradient;
                            frantic::volumetrics::levelset::trilerp_staggered_centered_signed_distance_gradient(
                                occlusionLS, faceCoord, occlusionGradient );

                            occlusionGradient.normalize();
                            vector3f occlusionNormal = -occlusionGradient;

                            vector3f fluidVelocity;
                            frantic::volumetrics::levelset::trilerp_vector3f_staggered(
                                velocityFieldRIS, inputStaggeredVelocityAcc, faceCoord, fluidVelocity );

                            vector3f relativeVelocity = fluidVelocity - occlusionVelocityAtFace;

                            if( vector3f::dot( relativeVelocity, occlusionNormal ) >= 0.f ) {
                                // we set the normal velocity component to the normal component of the occlusion
                                // velocity
                                vector3f normalComponent;
                                if( occlusionHasVelocity ) {
                                    normalComponent =
                                        vector3f::dot( occlusionVelocityAtFace, occlusionNormal ) * occlusionNormal;
                                }

                                vector3f tangentComponent =
                                    fluidVelocity - vector3f::dot( fluidVelocity, occlusionNormal ) * occlusionNormal;

                                resultVelocity[axis] = normalComponent[axis] + tangentComponent[axis];
                            } else {
                                // the fluid is unconstrained
                                resultVelocity[axis] = fluidVelocity[axis];
                            }
                        }
                    }

                } else {
                    // Will this face get a neumann boundary condition in the pressure solve ?
                    // If so, set it to the occlusion's velocity.
                    const bool faceInsideOcclusion =
                        occlusionLS[occlusionIndex] <= 0 || occlusionLS[adjOcclusionIndex] <= 0;
                    if( faceInsideOcclusion ) {
                        if( occlusionHasVelocity ) {
                            const float occlusionVelocityAtFace =
                                0.5f * ( occlusionVelocityAcc[occlusionIndex][axis] +
                                         occlusionVelocityAcc[adjOcclusionIndex][axis] );
                            resultVelocity[axis] = occlusionVelocityAtFace;
                        } else {
                            resultVelocity[axis] = 0;
                        }
                    }
                }
            }
        }

        outputStaggeredVelocity[velocityIndex] = resultVelocity;
    }

    memcpy( &inputStaggeredVelocityAcc[0], &outputStaggeredVelocity[0],
            inputStaggeredVelocityAcc.size() * sizeof( vector3f ) );
}

void frantic::fluids::apply_sim_boundary( const boundbox3& simBounds, rle_voxel_field& velocityField,
                                          BOUNDARY_CONDITION wallCondition, const frantic::tstring& faceStateChannel ) {
    // construct the face state field if required
    if( !velocityField.has_channel( faceStateChannel ) ) {
        throw std::runtime_error( "apply_sim_boundary() - Velocity Field does not have channel \"" +
                                  frantic::strings::to_string( faceStateChannel ) +
                                  "\" to hold the face boundary values" );
    }

    cout << "Provided Sim Bounds:" << simBounds << endl;

    rle_channel_accessor<boost::uint32_t> state =
        velocityField.get_channel_accessor<boost::uint32_t>( faceStateChannel );
    rle_channel_accessor<vector3f> velocity = velocityField.get_channel_accessor<vector3f>( _T("StaggeredVelocity") );

    FF_LOG( error ) << "**** (sim bounds )staggered velocity channel size=" << velocity.size() << " ****" << endl;

    boundary_policy_function_t boundaryPolicyFunction;
    switch( wallCondition ) {
    case NEUMANN:
        boundaryPolicyFunction = standard_solid_boundary_policy;
        break;
    case DIRICHLET:
        boundaryPolicyFunction = standard_free_boundary_policy;
        break;
    default:
        throw std::runtime_error( "apply_sim_boundary() - BOUNDARY_CONDITION " + lexical_cast<string>( wallCondition ) +
                                  " is not supported for the sim wall boundaries" );
    }

    const rle_index_spec& ris = velocityField.get_rle_index_spec();

    vector3 minCoord = simBounds.minimum();
    vector3 staggeredGhost = simBounds.maximum();
    vector3 maxCoord = staggeredGhost - vector3( 1 );

    if( boundaryPolicyFunction == standard_solid_boundary_policy ) {
        FF_LOG( debug ) << "Using standard SOLID boundary policy" << endl;
    }
    if( boundaryPolicyFunction == standard_free_boundary_policy ) {
        FF_LOG( debug ) << "Using standard FREE boundary policy" << endl;
    }

    boundbox3 voxelExtents( minCoord.x, staggeredGhost.x, minCoord.y, staggeredGhost.y, minCoord.z, minCoord.z );
    std::vector<boost::int32_t> dataIndices( voxelExtents.get_volume() );
    ris.fill_data_index_box( voxelExtents, &dataIndices[0] );

    int xsize = voxelExtents.xsize();
    int ysize = voxelExtents.ysize();
    int zsize = voxelExtents.zsize();
    for( int x = 0; x < xsize; x++ ) {
        for( int y = 0; y < ysize; ++y ) {
            int index = dataIndices[x + xsize * y];
            if( index >= 0 ) {
                if( y == ysize - 1 || x == xsize - 1 ) {
                    state[index] = 0;
                    velocity[index] = vector3f();
                    continue;
                }

                state[index] = set_face_boundary( state[index], 5, wallCondition );
                // TODO this should be done nicer, so that there is a clear way to set the velocity on the face
                // appropriately
                velocity[index].z = 0.f;
            }
        }
    }

    // the max z plane
    voxelExtents.set( minCoord.x, staggeredGhost.x, minCoord.y, staggeredGhost.y, maxCoord.z, maxCoord.z );
    ris.fill_data_index_box( voxelExtents, &dataIndices[0] );

    xsize = voxelExtents.xsize();
    ysize = voxelExtents.ysize();
    for( int x = 0; x < xsize; x++ ) {
        for( int y = 0; y < ysize; ++y ) {
            int index = dataIndices[x + xsize * y];
            if( index >= 0 ) {
                if( y == ysize - 1 || x == xsize - 1 ) {
                    state[index] = 0;
                    velocity[index] = vector3f();
                    continue;
                }

                state[index] = set_face_boundary( state[index], 4, wallCondition );
            }
        }
    }

    // the top ghost plane
    voxelExtents.set( minCoord.x, staggeredGhost.x, minCoord.y, staggeredGhost.y, staggeredGhost.z, staggeredGhost.z );
    ris.fill_data_index_box( voxelExtents, &dataIndices[0] );

    xsize = voxelExtents.xsize();
    ysize = voxelExtents.ysize();
    for( int x = 0; x < xsize; x++ ) {
        for( int y = 0; y < ysize; ++y ) {
            int index = dataIndices[x + xsize * y];
            if( index >= 0 ) {
                state[index] = 0;
                if( y == ysize - 1 || x == xsize - 1 ) {
                    velocity[index] = vector3f();
                    continue;
                }

                state[index] = set_face_boundary( state[index], 5, wallCondition );
                if( wallCondition == NEUMANN )
                    velocity[index].z = 0.f;
            }
        }
    }

    voxelExtents.set( minCoord.x, staggeredGhost.x, minCoord.y, minCoord.y, minCoord.z, staggeredGhost.z );
    dataIndices.resize( voxelExtents.get_volume() );
    ris.fill_data_index_box( voxelExtents, &dataIndices[0] );

    xsize = voxelExtents.xsize();
    ysize = voxelExtents.ysize();
    zsize = voxelExtents.zsize();

    for( int z = 0; z < zsize; ++z ) {
        for( int y = 0; y < ysize; ++y ) {
            for( int x = 0; x < xsize; x++ ) {

                int index = dataIndices[x + xsize * ( y + ysize * z )];
                if( index >= 0 ) {
                    if( z == zsize - 1 || x == xsize - 1 ) {
                        state[index] = set_face_boundary( state[index], 3, INVALID );
                        velocity[index] = vector3f();
                    } else {
                        state[index] = set_face_boundary( state[index], 3, wallCondition );
                        // TODO this should be done nicer, so that there is a clear way to set the velocity on the face
                        // appropriately
                        velocity[index].y = 0.f;
                    }
                }
            }
        }
    }

    voxelExtents.set( minCoord.x, staggeredGhost.x, maxCoord.y, maxCoord.y, minCoord.z, staggeredGhost.z );
    dataIndices.resize( voxelExtents.get_volume() );
    ris.fill_data_index_box( voxelExtents, &dataIndices[0] );

    xsize = voxelExtents.xsize();
    ysize = voxelExtents.ysize();
    zsize = voxelExtents.zsize();

    for( int z = 0; z < zsize; ++z ) {
        for( int y = 0; y < ysize; ++y ) {
            for( int x = 0; x < xsize; x++ ) {

                int index = dataIndices[x + xsize * ( y + ysize * z )];
                if( index >= 0 ) {
                    if( z == zsize - 1 || x == xsize - 1 ) {
                        state[index] = set_face_boundary( state[index], 2, INVALID );
                        velocity[index] = vector3f();
                    } else {
                        state[index] = set_face_boundary( state[index], 2, wallCondition );
                    }
                }
            }
        }
    }

    voxelExtents.set( minCoord.x, staggeredGhost.x, staggeredGhost.y, staggeredGhost.y, minCoord.z, staggeredGhost.z );
    dataIndices.resize( voxelExtents.get_volume() );
    ris.fill_data_index_box( voxelExtents, &dataIndices[0] );

    xsize = voxelExtents.xsize();
    ysize = voxelExtents.ysize();
    zsize = voxelExtents.zsize();

    for( int z = 0; z < zsize; ++z ) {
        for( int y = 0; y < ysize; ++y ) {
            for( int x = 0; x < xsize; x++ ) {

                int index = dataIndices[x + xsize * ( y + ysize * z )];
                if( index >= 0 ) {
                    state[index] = 0;
                    if( z == zsize - 1 || x == xsize - 1 ) {
                        state[index] = set_face_boundary( state[index], 3, INVALID );
                        velocity[index] = vector3f();
                    } else {
                        state[index] = set_face_boundary( state[index], 3, wallCondition );
                        if( wallCondition == NEUMANN )
                            velocity[index].y = 0.f;
                    }
                }
            }
        }
    }

    voxelExtents.set( minCoord.x, minCoord.x, minCoord.y, staggeredGhost.y, minCoord.z, staggeredGhost.z );
    dataIndices.resize( voxelExtents.get_volume() );
    ris.fill_data_index_box( voxelExtents, &dataIndices[0] );

    xsize = voxelExtents.xsize();
    ysize = voxelExtents.ysize();
    zsize = voxelExtents.zsize();

    for( int z = 0; z < zsize; ++z ) {
        for( int y = 0; y < ysize; ++y ) {
            for( int x = 0; x < xsize; x++ ) {

                int index = dataIndices[x + xsize * ( y + ysize * z )];
                if( index >= 0 ) {
                    if( z == zsize - 1 || y == ysize - 1 ) {
                        state[index] = set_face_boundary( state[index], 1, INVALID );
                        velocity[index] = vector3f();

                    } else {
                        state[index] = set_face_boundary( state[index], 1, wallCondition );
                        velocity[index].x = 0.f;
                    }
                }
            }
        }
    }

    voxelExtents.set( maxCoord.x, maxCoord.x, minCoord.y, staggeredGhost.y, minCoord.z, staggeredGhost.z );
    dataIndices.resize( voxelExtents.get_volume() );
    ris.fill_data_index_box( voxelExtents, &dataIndices[0] );

    xsize = voxelExtents.xsize();
    ysize = voxelExtents.ysize();
    zsize = voxelExtents.zsize();

    for( int z = 0; z < zsize; ++z ) {
        for( int y = 0; y < ysize; ++y ) {
            for( int x = 0; x < xsize; x++ ) {

                int index = dataIndices[x + xsize * ( y + ysize * z )];
                if( index >= 0 ) {
                    if( z == zsize - 1 || y == ysize - 1 ) {
                        velocity[index] = vector3f();
                        state[index] = set_face_boundary( state[index], 0, INVALID );
                    } else {
                        state[index] = set_face_boundary( state[index], 0, wallCondition );
                    }
                }
            }
        }
    }

    voxelExtents.set( staggeredGhost.x, staggeredGhost.x, minCoord.y, staggeredGhost.y, minCoord.z, staggeredGhost.z );
    dataIndices.resize( voxelExtents.get_volume() );
    ris.fill_data_index_box( voxelExtents, &dataIndices[0] );

    xsize = voxelExtents.xsize();
    ysize = voxelExtents.ysize();
    zsize = voxelExtents.zsize();

    for( int z = 0; z < zsize; ++z ) {
        for( int y = 0; y < ysize; ++y ) {
            for( int x = 0; x < xsize; x++ ) {

                int index = dataIndices[x + xsize * ( y + ysize * z )];
                if( index >= 0 ) {
                    state[index] = 0;

                    if( z == zsize - 1 || y == ysize - 1 ) {
                        state[index] = set_face_boundary( state[index], 1, INVALID );
                        velocity[index] = vector3f();
                    } else {
                        state[index] = set_face_boundary( state[index], 1, wallCondition );
                        if( wallCondition == NEUMANN )
                            velocity[index].x = 0.f;
                    }
                }
            }
        }
    }
}

void frantic::fluids::zero_sim_boundary_velocity( const boundbox3& voxelBounds, rle_voxel_field& velocityField,
                                                  const frantic::tstring& staggeredVelocityChannelName,
                                                  const frantic::tstring& staggeredPopulatedChannelName ) {
    if( velocityField.size() == 0 )
        return;

    if( !velocityField.has_channel( staggeredVelocityChannelName ) ) {
        throw std::runtime_error( "zero_sim_boundary_velocity() Error: the voxel field is missing the specified "
                                  "staggered velocity channel \'" +
                                  frantic::strings::to_string( staggeredVelocityChannelName ) + "\'" );
    }
    rle_channel_accessor<vector3f> staggeredVelocityAcc =
        velocityField.get_channel_accessor<vector3f>( staggeredVelocityChannelName );

    boost::uint8_t* staggeredPopulatedChannel = NULL;
    if( !staggeredPopulatedChannelName.empty() ) {
        if( !velocityField.has_channel( staggeredPopulatedChannelName ) ) {
            throw std::runtime_error( "zero_sim_boundary_velocity() Error: the voxel field is missing the specified "
                                      "staggered populated channel \'" +
                                      frantic::strings::to_string( staggeredPopulatedChannelName ) + "\'" );
        } else {
            rle_channel_accessor<boost::uint8_t> staggeredPopulatedAcc =
                velocityField.get_channel_accessor<unsigned char>( staggeredPopulatedChannelName );
            if( !staggeredPopulatedAcc.size() ) {
                throw std::runtime_error(
                    "zero_sim_boundary_velocity() Error: the velocity field has size > 0, but the "
                    "staggered populated channel has size 0." );
            }
            staggeredPopulatedChannel = &staggeredPopulatedAcc[0];
        }
    }

    const rle_index_spec& ris = velocityField.get_rle_index_spec();

    const vector3 minCoord = voxelBounds.minimum();
    const vector3 staggeredGhost = voxelBounds.maximum();

    // the staggered populated bit for each axis
    const boost::uint8_t staggeredPopulatedBit[] = { 0x01, 0x02, 0x04 };
    // const boost::uint8_t staggeredPopulatedAll = 0x07;

    std::vector<boost::int32_t> dataIndices;
    boundbox3 voxelExtents;

    for( int face = 0; face < 6; ++face ) {
        switch( face ) {
        case rae_index_x_pos:
            voxelExtents.set( staggeredGhost.x, staggeredGhost.x, minCoord.y, staggeredGhost.y, minCoord.z,
                              staggeredGhost.z );
            break;
        case rae_index_x_neg:
            voxelExtents.set( minCoord.x, minCoord.x, minCoord.y, staggeredGhost.y, minCoord.z, staggeredGhost.z );
            break;
        case rae_index_y_pos:
            voxelExtents.set( minCoord.x, staggeredGhost.x, staggeredGhost.y, staggeredGhost.y, minCoord.z,
                              staggeredGhost.z );
            break;
        case rae_index_y_neg:
            voxelExtents.set( minCoord.x, staggeredGhost.x, minCoord.y, minCoord.y, minCoord.z, staggeredGhost.z );
            break;
        case rae_index_z_pos:
            voxelExtents.set( minCoord.x, staggeredGhost.x, minCoord.y, staggeredGhost.y, staggeredGhost.z,
                              staggeredGhost.z );
            break;
        case rae_index_z_neg:
            voxelExtents.set( minCoord.x, staggeredGhost.x, minCoord.y, staggeredGhost.y, minCoord.z, minCoord.z );
            break;
        default:
            throw std::runtime_error( "zero_sim_boundary_velocity() Internal Error: invalid face number " +
                                      boost::lexical_cast<std::string>( face ) );
        }

        const int axis = neighbor_index_axis( face );
        // const bool isGhostVoxel = is_neighbor_index_direction_positive( face );

        const int indexCount = voxelExtents.get_volume();

        dataIndices.resize( indexCount );
        ris.fill_data_index_box( voxelExtents, &dataIndices[0] );

        BOOST_FOREACH( const boost::int32_t dataIndex, dataIndices ) {
            if( dataIndex >= 0 ) {
                staggeredVelocityAcc[dataIndex][axis] = 0;
                if( staggeredPopulatedChannel ) {
                    staggeredPopulatedChannel[dataIndex] |= staggeredPopulatedBit[axis];
                }
            }
        }
    }
}

bool frantic::fluids::standard_free_boundary_policy( const boundbox3& bounds, const vector3& voxelCoord, int face,
                                                     boost::uint32_t& outValue ) {
    // the sim bounds overrides a fluid condition but the occlusion condition overrides the free boundary
    // condition so the occlusion test should be last
    bool set = false;
    if( bounds.contains( voxelCoord ) ) {
        if( face % 2 != 0 ) {
            if( bounds.minimum()[face / 2] == voxelCoord[face / 2] ) {
                // FF_LOG(debug) << "\t\tboundary cell: min" << endl;
                outValue = set_face_boundary( outValue, face, DIRICHLET );
                set = true;
            }
            // check that the coordinate is not in a 'ghost' cell (which stores the outside face data)
            if( bounds.maximum()[face / 2] == voxelCoord[face / 2] ) {
                // FF_LOG(debug) << "\t\tboundary cell: ghost" << endl;
                outValue = set_face_boundary( outValue, face, DIRICHLET );
                set = true;
            }
        } else {
            vector3 realMax = bounds.maximum() - vector3( 1 );

            if( realMax[face / 2] == voxelCoord[face / 2] ) {
                // FF_LOG(debug) << "\t\tboundary cell: max" << endl;
                outValue = set_face_boundary( outValue, face, DIRICHLET );
                set = true;
            }
        }
    }
    return set;
}

bool frantic::fluids::standard_solid_boundary_policy( const boundbox3& bounds, const vector3& voxelCoord, int face,
                                                      boost::uint32_t& outValue ) {
    // the sim bounds overrides a fluid condition but the occlusion condition overrides the free boundary
    // condition so the occlusion test should be last
    frantic::logging::set_logging_level( 5 );
    bool set = false;
    if( face % 2 != 0 ) { // its a neg face
        // cout << "boundary testing face " << face << " voxel " <<  voxelCoord << endl;
        if( bounds.minimum()[face / 2] == voxelCoord[face / 2] ) {
            FF_LOG( debug ) << "\t\tboundary cell: min" << endl;
            outValue = set_face_boundary( outValue, face, NEUMANN );
            set = true;
        }
        // check that the coordinate is not in a 'ghost' cell (which stores the outside face data)
        if( bounds.maximum()[face / 2] == voxelCoord[face / 2] ) {
            FF_LOG( debug ) << "\t\tboundary cell: ghost" << endl;
            outValue = set_face_boundary( outValue, face, NEUMANN );
            //// and set all other faces to invalid for safety's sake
            // for( int i =0; i<6; ++i){
            //	if( i!=face )
            //		outValue = set_face_boundary(outValue,i,INVALID);
            // }

            set = true;
        }
    } else {
        // cout << "boundary testing face " << face << " voxel " <<  voxelCoord << endl;
        //  its a pos face

        // check a positive face against the real max (to skip the ghost cell)
        vector3 realMax = bounds.maximum() - vector3( 1 );

        if( realMax[face / 2] == voxelCoord[face / 2] ) {
            FF_LOG( debug ) << "\t\tboundary cell: max" << endl;
            outValue = set_face_boundary( outValue, face, NEUMANN );
            set = true;
        }

        // we want to know if the next voxel over is the boundary so we can correctly close it off
        if( bounds.minimum()[face / 2] == voxelCoord[face / 2] + 1 ) {
            outValue = set_face_boundary( outValue, face, NEUMANN );
            set = true;
        }
    }

    frantic::logging::set_logging_level( 3 );
    return set;
}
