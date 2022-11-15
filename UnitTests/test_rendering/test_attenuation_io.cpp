// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/assign.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/join.hpp>

#include <frantic/files/scoped_file_cleanup.hpp>

#include <frantic/rendering/deep_attenuation_loaders.hpp>
#include <frantic/rendering/deep_attenuation_savers.hpp>

using namespace frantic::graphics;

TEST( attenuation_io, cubeface_bilinear_sampling ) {

    static const std::vector<vector3f> negLocations = boost::assign::list_of( vector3f( -1.f, 0.f, 0.f ) )(
        vector3f( 0.f, -1.f, 0.f ) )( vector3f( 0.f, 0.f, -1.f ) );

    static const std::vector<vector3f> posLocations =
        boost::assign::list_of( vector3f( 1.f, 0.f, 0.f ) )( vector3f( 0.f, 1.f, 0.f ) )( vector3f( 0.f, 0.f, 1.f ) );

    static const boost::range::joined_range<const std::vector<vector3f>, const std::vector<vector3f>> rng =
        boost::join( negLocations, posLocations );
    static const std::vector<vector3f> allLocations( rng.begin(), rng.end() );

    static const std::vector<std::vector<vector3f>> locationSets =
        boost::assign::list_of( allLocations )( negLocations )( posLocations );

    static const std::vector<std::vector<vector3f>> excludedSets =
        boost::assign::list_of( std::vector<vector3f>() )( posLocations )( negLocations );

    boost::filesystem::path map_path =
        boost::filesystem::unique_path( boost::filesystem::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.exr" );
    frantic::files::scoped_file_cleanup fileCleanup;
    fileCleanup.add( map_path );

    int i = 0;

    for( std::vector<std::vector<vector3f>>::const_iterator setIt = locationSets.begin(), exIt = excludedSets.begin();
         setIt != locationSets.end(); ++setIt, ++exIt, ++i ) {
        std::cout << "Iteration: " << i << std::endl;

        {
            frantic::rendering::cubeface_atten_exr_saver saver( map_path.string<frantic::tstring>(), 512, 1, 1.f,
                                                                false );
            for( std::vector<vector3f>::const_iterator it = setIt->begin(); it != setIt->end(); ++it ) {
                saver.add_sample_bilinear( *it, alpha3f( 1.f ) );
            }
            saver.write_file();
        }

        {
            frantic::rendering::cubeface_atten_exr_loader loader( map_path.string<frantic::tstring>() );
            for( std::vector<vector3f>::const_iterator it = setIt->begin(); it != setIt->end(); ++it ) {
                std::cout << *it << std::endl;
                EXPECT_EQ( 1.f, loader.get_zdepth_bilinear( *it ) );
                EXPECT_EQ( alpha3f( 1.f ), loader.get_sample_bilinear( 2.f * *it ) );
            }

            for( std::vector<vector3f>::const_iterator it = exIt->begin(); it != exIt->end(); ++it ) {
                EXPECT_EQ( std::numeric_limits<float>::max(), loader.get_zdepth_bilinear( *it ) );
                EXPECT_EQ( alpha3f( 0.f ), loader.get_sample_bilinear( 2.f * *it ) );
            }
        }
    }
}
