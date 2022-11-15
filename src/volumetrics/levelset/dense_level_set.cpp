// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/volumetrics/levelset/dense_level_set.hpp>

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::channels;
using namespace frantic::graphics;
using namespace frantic::graphics2d;
using namespace frantic::volumetrics;

void frantic::volumetrics::levelset::dense_level_set::dump( std::ostream& out ) const {
    out << "DUMPING DENSE LEVEL SET\n";
    out << "Voxel count in bounds: " << m_voxelBoundsSize.volume() << "\n";
    out << "Outer bounds: " << outer_bounds() << "\n";

    out << "Voxel Coordinate System: " << m_voxelCoordSystem << "\n";
    out << "\n";
    out << "Number of named channels: " << m_namedChannels.size() << "\n";
    int index = 0;
    for( map<frantic::tstring, dense_level_set_channel>::const_iterator i = m_namedChannels.begin();
         i != m_namedChannels.end(); ++i ) {
        out << "Channel " << index++ << ", \"" << frantic::strings::to_string( i->first ) << "\"\n";
        out << "Data type: " << i->second.type_str() << "\n";
        // out << "Data: ";
        // channels::channel_data_type_print( out, ",", i->second.arity() * i->second.size(), i->second.data_type(),
        // i->second.data() ); out << "\n";
    }
    out << "FINISHED DUMPING DENSE LEVEL SET" << endl;
}

void frantic::volumetrics::levelset::dense_level_set::fill_plane( const boundrect2& voxelXYExtents, int voxelZ,
                                                                  std::vector<float>& outVoxelCornerValues ) const {
    fill_box( boundbox3( vector3( voxelXYExtents.minimum().x, voxelXYExtents.minimum().y, voxelZ ),
                         vector3( voxelXYExtents.maximum().x, voxelXYExtents.maximum().y, voxelZ ) ),
              outVoxelCornerValues );
}

// This function has only been tested as used through fill_plane so far.  Need to validate that it works in more general
// settings.
void frantic::volumetrics::levelset::dense_level_set::fill_box( const frantic::graphics::boundbox3& /*voxelExtents*/,
                                                                std::vector<float>& /*outVoxelCornerValues*/ ) const {
    throw runtime_error( "frantic::volumetrics::levelset::dense_level_set::fill_box unimplemented" );
    // TODO
}

// This sets the rle level set to the given data, using swap functions to set the data efficiently.  Note that
// the variables you pass in will have arbitrary values in them after this function is called.
void frantic::volumetrics::levelset::dense_level_set::set_with_swap( voxel_coord_system& vcs,
                                                                     const graphics::boundbox3& voxelBounds,
                                                                     std::vector<float>& distanceData ) {
    if( voxelBounds.get_volume() != (int)distanceData.size() )
        throw runtime_error( "dense_level_set.set_with_swap: The size of the data provided, " +
                             boost::lexical_cast<std::string>( distanceData.size() ) +
                             ", doesn't match the data size of the RLE Index Spec provided, " +
                             boost::lexical_cast<std::string>( voxelBounds.get_volume() ) + "." );

    m_voxelCoordSystem.swap( vcs );
    m_voxelBoundsMin = voxelBounds.minimum();
    m_voxelBoundsSize = voxelBounds.size();
    m_distanceData.swap( distanceData );
    m_namedChannels.clear();
}
