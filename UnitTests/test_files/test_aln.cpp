// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <boost/foreach.hpp>
#include <frantic/files/aln.hpp>

#include "gtest-helper.h"

using namespace frantic;
using namespace frantic::files;
using frantic::tstring;
using frantic::graphics::transform4f;
using frantic::graphics::transform4fd;
using frantic::graphics::vector3fd;
using std::make_pair;

TEST( ALN, ParseFileSimple ) {
    std::vector<std::pair<tstring, transform4fd>> testData;
    testData.push_back(
        make_pair( _T("hdri_colors.csv"), transform4fd( -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 ) ) );
    testData.push_back(
        make_pair( _T("singleplusoffset.csv"), transform4fd( -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 2, -4, 7, 1 ) ) );

    {
        std::vector<std::pair<tstring, transform4fd>> actualTransforms;
        aln::get_transforms<transform4fd>( _T("TestInputs/align.aln"), actualTransforms );
        ASSERT_EQ( testData.size(), actualTransforms.size() );
        for( int i = 0; i < testData.size(); ++i ) {
            EXPECT_EQ( testData[i].first, actualTransforms[i].first );
            EXPECT_TRANSFORM4FD_EQ( testData[i].second, actualTransforms[i].second );
        }
    }
    {
        std::vector<std::pair<tstring, transform4f>> actualTransforms;
        aln::get_transforms<transform4f>( _T("TestInputs/align.aln"), actualTransforms );
        ASSERT_EQ( testData.size(), actualTransforms.size() );
        for( int i = 0; i < testData.size(); ++i ) {
            EXPECT_EQ( testData[i].first, actualTransforms[i].first );
            EXPECT_TRANSFORM4F_EQ( transform4f( testData[i].second ), actualTransforms[i].second );
        }
    }
}

#if defined( LZ4_AVAILABLE ) && defined( E57_AVAILABLE )
TEST( ALN, ParseFileVariety ) {
    std::vector<std::pair<tstring, transform4fd>> testData;
    testData.push_back(
        make_pair( _T("intColor.prt"),
                   transform4fd( 717.1969444082445, 326.69458919274257, -647.617699401002, -406.9242697290574,
                                 954.5372303016541, 209.8710170122074, -930.1021600843956, 901.5319708453505,
                                 -121.32680845462392, 490.7905621404634, -585.677279722057, -210.86889162471368,
                                 501.4218330272854, 272.31519910961765, 859.2444155562719, 116.40924282311926 ) ) );
    testData.push_back(
        make_pair( _T("scale_offset_las.las"),
                   transform4fd( -638.4616140774733, -43.983628320985304, 501.4727786741853, 486.53724482047346,
                                 -775.6111858930155, 224.62615588950052, 324.721600726783, -929.2668937895945,
                                 -433.32686291981975, 358.8467232765372, 764.658555590504, -300.26329198458154,
                                 868.0703744984751, 231.61573297133532, 296.9875966928348, -686.9035639865401 ) ) );
    testData.push_back(
        make_pair( _T("pts_inconsistent_channels.pts"),
                   transform4fd( -70.50402526766368, -96.55087768653425, -428.37609341833956, -220.5565932560669,
                                 216.03507836695394, -386.89903398933984, 373.9322301850568, 625.5427621271613,
                                 -394.20353572310887, -22.323114007633194, 11.17992561891333, -358.93937542514334,
                                 -235.09740818014245, 286.7448485441296, 845.6884294431557, -508.99853381760596 ) ) );
    testData.push_back(
        make_pair( _T("bunnyFloat_wlod.sprt"),
                   transform4fd( -143.9299971833252, -900.0755749296745, 88.8513525753383, 924.2281160751047,
                                 463.4312547106458, -819.9498082622945, -207.56634201010309, -435.8655640213933,
                                 363.49353131727753, -38.055306469356765, 979.0353561480686, 55.527618719156635,
                                 -749.8229199248912, 226.70452494176152, 103.71354414726375, 154.06840964133153 ) ) );
    testData.push_back(
        make_pair( _T("hdri_colors.csv"),
                   transform4fd( -937.114160904982, 907.8036620211371, -661.0550732641842, 5.453083443112973,
                                 -810.5345416466141, 444.69723455904114, 469.54638840798316, -637.8326506154287,
                                 -76.66969015601398, -635.1160898894836, 237.93248665477563, 282.4855555695624,
                                 74.52991470102756, -317.3180783159828, -937.144805229146, 809.4780975471595 ) ) );
    testData.push_back(
        make_pair( _T("wnormals.e57"),
                   transform4fd( -891.1280813532369, 758.9904178188804, 93.53572629705718, -745.4208916778937,
                                 -181.6040461628603, 923.9285520155013, 787.7624143398384, -223.2669876565734,
                                 -337.0265665480101, 504.08057930309224, -697.6268496448448, 221.82619107676896,
                                 500.2995113067045, 612.9204821600688, -76.94768439802704, 532.0635043106622 ) ) );
    for( int i = 0; i < testData.size(); ++i )
        testData[i].second.transpose();

    {
        std::vector<std::pair<tstring, transform4fd>> actualTransforms;
        aln::get_transforms<transform4fd>( _T("TestInputs/variety.aln"), actualTransforms );
        ASSERT_EQ( testData.size(), actualTransforms.size() );
        for( int i = 0; i < testData.size(); ++i ) {
            EXPECT_EQ( testData[i].first, actualTransforms[i].first );
            EXPECT_TRANSFORM4FD_EQ( testData[i].second, actualTransforms[i].second );
        }
    }
    {
        std::vector<std::pair<tstring, transform4f>> actualTransforms;
        aln::get_transforms<transform4f>( _T("TestInputs/variety.aln"), actualTransforms );
        ASSERT_EQ( testData.size(), actualTransforms.size() );
        for( int i = 0; i < testData.size(); ++i ) {
            EXPECT_EQ( testData[i].first, actualTransforms[i].first );
            EXPECT_TRANSFORM4F_EQ( transform4f( testData[i].second ), actualTransforms[i].second );
        }
    }
}
#endif

TEST( ALN, TooManyScans ) {
    typedef std::vector<std::pair<tstring, transform4fd>> mapping_t;
    ASSERT_THROW(
        {
            mapping_t actualTransforms;
            aln::get_transforms<transform4fd>( _T("TestInputs/aln_toomany.aln"), actualTransforms );
        },
        std::runtime_error );
}

TEST( ALN, TooFewScans ) {
    typedef std::vector<std::pair<tstring, transform4fd>> mapping_t;
    ASSERT_THROW(
        {
            mapping_t actualTransforms;
            aln::get_transforms<transform4fd>( _T("TestInputs/aln_toofew.aln"), actualTransforms );
        },
        std::runtime_error );
}

TEST( ALN, TransformScan ) {
    vector3fd expectedA[] = { // hdri_colors.csv, reflected across x
                              vector3fd( 0.0, 0.0, 0.0 ), vector3fd( -1.0, 0.0, 0.0 ), vector3fd( 0.0, 1.0, 0.0 ),
                              vector3fd( 0.0, 0.0, 1.0 ),
                              // singleplusoffset.csv, reflected across (x, y, z), then translated (2, -4, 7)
                              vector3fd( -100000000.0 + 2, -400.0 - 4, -401000001.0 + 7 ),
                              vector3fd( -100000001.0 + 2, -500.0 - 4, -401000000.0 + 7 ),
                              vector3fd( -100000010.0 + 2, -600.0 - 4, -401000002.0 + 7 ),
                              vector3fd( -100000100.0 + 2, -700.0 - 4, -401000003.0 + 7 ),
                              vector3fd( -100000000.0 + 2, -900.0 - 4, -401000004.0 + 7 ) };
    std::vector<vector3fd> expectedPositions;
    expectedPositions.assign( expectedA, expectedA + 9 );

    particles::particle_file_stream_factory_object factory;
    std::vector<vector3fd> actualPositions;
    factory.set_position_type_hint( channels::data_type_float64 );
    particles::streams::particle_istream_ptr is = factory.create_istream( _T("TestInputs/align.aln") );

    particles::particle_array pArr( is->get_channel_map() );
    pArr.insert_particles( is );

    const channels::channel_const_cvt_accessor<vector3fd> accessor =
        pArr.get_channel_map().get_const_cvt_accessor<vector3fd>( _T("Position") );
    for( int i = 0; i < pArr.size(); ++i )
        actualPositions.push_back( accessor( pArr[i] ) );

    ASSERT_EQ( expectedPositions.size(), actualPositions.size() );
    for( int i = 0; i < expectedPositions.size(); ++i ) {
        EXPECT_VECTOR3FD_EQ( expectedPositions[i], actualPositions[i] );
    }
}
