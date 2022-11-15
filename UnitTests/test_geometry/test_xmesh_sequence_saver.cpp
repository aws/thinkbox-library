// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/files/files.hpp>
#include <frantic/files/scoped_file_cleanup.hpp>
#include <frantic/geometry/xmesh_sequence_saver.hpp>
#include <frantic/locale/locale.hpp>

class XMeshSequenceSaverTestFixture : public ::testing::Test {
  public:
    virtual void SetUp() {
        using namespace boost::filesystem;

        m_tempDir = temp_directory_path() / unique_path();
        create_directory( m_tempDir );
        m_cleanup.add( m_tempDir );
    }

    const boost::filesystem::path& get_temp_dir() const { return m_tempDir; }

  private:
    boost::filesystem::path m_tempDir;
    frantic::files::scoped_file_cleanup m_cleanup;
};

TEST_F( XMeshSequenceSaverTestFixture, SpanishLocale ) {
    using namespace frantic::files;
    using namespace frantic::geometry;
    using namespace frantic::graphics;
    using namespace frantic::locale;

    const double lengthUnitScale = 0.1;

    const frantic::tstring filename = to_tstring( get_temp_dir() / _T("Spanish.xmesh") );

    { // scope for Spanish locale
#ifdef _WIN32
#if defined( _MSC_VER ) && _MSC_VER >= 1700
        const char* spanish = "es-ES";
#else
        const char* spanish = "Spanish";
#endif
#else
        const char* spanish = "es_ES.UTF-8";
#endif
        set_locale_in_scope setLocale( spanish );

        { // scope for writing
            xmesh_metadata metadata;
            metadata.set_length_unit( lengthUnitScale, xmesh_metadata::length_unit_meters );

            trimesh3 mesh;
            mesh.set_to_box( boundbox3f( vector3f( 0 ), vector3f( 0.1f ) ) );

            frantic::geometry::xmesh_sequence_saver xss;
            xss.write_xmesh( mesh, metadata, filename );
        }

        { // scope for reading in Spanish locale
            xmesh_metadata metadata;
            read_xmesh_metadata( filename, metadata );
            ASSERT_TRUE( metadata.has_boundbox() );
            EXPECT_FLOAT_EQ( 0, metadata.get_boundbox().xminimum() );
            EXPECT_FLOAT_EQ( 0, metadata.get_boundbox().yminimum() );
            EXPECT_FLOAT_EQ( 0, metadata.get_boundbox().zminimum() );
            EXPECT_FLOAT_EQ( 0.1f, metadata.get_boundbox().xmaximum() );
            EXPECT_FLOAT_EQ( 0.1f, metadata.get_boundbox().ymaximum() );
            EXPECT_FLOAT_EQ( 0.1f, metadata.get_boundbox().zmaximum() );

            EXPECT_DOUBLE_EQ( lengthUnitScale, metadata.get_length_unit_scale() );
        }
    }

    { // scope for reading in "C" locale
        set_locale_in_scope locale( "C" );

        xmesh_metadata metadata;
        read_xmesh_metadata( filename, metadata );
        ASSERT_TRUE( metadata.has_boundbox() );
        EXPECT_FLOAT_EQ( 0, metadata.get_boundbox().xminimum() );
        EXPECT_FLOAT_EQ( 0, metadata.get_boundbox().yminimum() );
        EXPECT_FLOAT_EQ( 0, metadata.get_boundbox().zminimum() );
        EXPECT_FLOAT_EQ( 0.1f, metadata.get_boundbox().xmaximum() );
        EXPECT_FLOAT_EQ( 0.1f, metadata.get_boundbox().ymaximum() );
        EXPECT_FLOAT_EQ( 0.1f, metadata.get_boundbox().zmaximum() );

        EXPECT_DOUBLE_EQ( lengthUnitScale, metadata.get_length_unit_scale() );
    }
}
