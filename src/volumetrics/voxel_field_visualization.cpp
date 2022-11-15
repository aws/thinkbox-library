// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on
#include <frantic/volumetrics/voxel_field_visualization.hpp>

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::graphics;
using namespace frantic::geometry;
using namespace frantic::volumetrics::levelset;

void frantic::volumetrics::visualization::convert_voxel_channel_to_debug_mesh(
    const frantic::volumetrics::voxel_coord_system& vcs, const frantic::volumetrics::levelset::rle_index_spec& ris,
    float cubeSize, const_rle_channel_general_accessor& channelAcc,
    frantic::volumetrics::visualization::debug_mesh_color_policy& colorPolicy, frantic::geometry::trimesh3& outMesh ) {
    trimesh3 result, tempMesh;
    frantic::tstring colorChannelName = _T("Color");

    // Initialize the tempMesh to a box.  For performance reasons, this box is never reallocated, it just has its
    // vertices moved about.
    tempMesh.set_to_exploded_box();

    // Add a color channel, and get an accessor to it
    tempMesh.add_vertex_channel<color3f>( colorChannelName );
    trimesh3_vertex_channel_accessor<color3f> ca = tempMesh.get_vertex_channel_accessor<color3f>( colorChannelName );

    // vector<string> names;
    // field.get_channel_names( names );
    // cout<< "field has channels: ";
    // for(int i=0; i!=names.size(); ++i){
    //	cout << names[i] << endl;
    // }

    // const_rle_channel_general_accessor channelAcc = field.get_channel_general_accessor( stateChannelName );

    // const_rle_channel_accessor<boost::uint32_t> state = field.get_channel_accessor<boost::uint32_t>("FaceState");

    // const voxel_coord_system& vcs = field.get_voxel_coord_system();
    // const rle_index_spec& ris = field.get_rle_index_spec();

    for( rle_defined_and_adj_iterator i( ris, false ), ie( ris, true ); i != ie; ++i ) {
        boundbox3f bounds( vcs.get_world_voxel_center( i.get_coord() ) );
        bounds.expand( 0.5f * cubeSize * vcs.voxel_length() );

        tempMesh.set_existing_exploded_box_verts( bounds );
        color3f col( 0 );

        // vector<color3f> colors(24);

        // if the policy says to add the cube
        if( colorPolicy.compute_vertex_colors( i, channelAcc, &ca[0] ) ) {

            // add the cube
            result.combine( tempMesh );
        }

        /*	boost::uint32_t faceValues = state[i.get_data_index()];

          int reorderedFaces[6] = {1,0,3,2,5,4};
          int invalidCount = 0;

          for( int face=0, v=0;face<6;++face) {
            boost::uint32_t value = (faceValues&frantic::fluids::face_masks::get_mask(reorderedFaces[face])) >>
          frantic::fluids::face_masks::get_shift(reorderedFaces[face]); switch( value  ) { case
          frantic::fluids::INVALID: col = color3f(0.f,0.8f,1.f);
                ++invalidCount;
                break;
              case frantic::fluids::NONE:
                col = color3f(0.f,0.f,1.f);
                break;
              case frantic::fluids::DIRICHLET:
                col = color3f(0.f,1.f,0.f);
                break;
              case frantic::fluids::NEUMANN:
                col = color3f(1.f,0.f,0.f);
                break;
              default:
                throw std::runtime_error("For voxel " + i.get_coord().str() + " the boundary condition value for face "
          + lexical_cast<std::string>(face) + " is not one of the set values (0,1,2,4) it is " +
          lexical_cast<string>(value)
          );
            }

            for( size_t i = v; i != v+4; ++i ) {
              ca[i] = col;
            }
            v+=4;
          }

          if( invalidCount != 6 )
            result.combine(tempMesh);*/
    }

    result.swap( outMesh );
}
