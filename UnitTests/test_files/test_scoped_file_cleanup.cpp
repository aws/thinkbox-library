// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/files/scoped_file_cleanup.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include <fstream>

TEST( ScopedFileCleanup, SimpleSingleFile ) {
    using namespace frantic::files;

    boost::filesystem::path testPath = boost::filesystem::unique_path(
        boost::filesystem::temp_directory_path() / boost::filesystem::path( "%%%%-%%%%-%%%%-%%%%.tmp" ) );

    EXPECT_FALSE( boost::filesystem::exists( testPath ) );

    {
        scoped_file_cleanup fileCleanup;
        fileCleanup.add( testPath );

        std::ofstream os( testPath.c_str() );

        os << "Lorem ipsum";

        os.close();

        EXPECT_TRUE( boost::filesystem::exists( testPath ) );
    }

    EXPECT_FALSE( boost::filesystem::exists( testPath ) );
}

TEST( ScopedFileCleanup, SimpleNonExistantFile ) {
    using namespace frantic::files;

    boost::filesystem::path testPath = boost::filesystem::unique_path(
        boost::filesystem::temp_directory_path().wstring() / boost::filesystem::path( "%%%%-%%%%-%%%%-%%%%.tmp" ) );

    EXPECT_FALSE( boost::filesystem::exists( testPath ) );

    {
        scoped_file_cleanup fileCleanup;
        fileCleanup.add( testPath );
    }

    EXPECT_FALSE( boost::filesystem::exists( testPath ) );
}

TEST( ScopedFileCleanup, CleanEmptyDirectory ) {
    using namespace frantic::files;

    boost::filesystem::path testPath = boost::filesystem::unique_path(
        boost::filesystem::temp_directory_path().wstring() / boost::filesystem::path( "%%%%-%%%%-%%%%-%%%%" ) );

    EXPECT_FALSE( boost::filesystem::exists( testPath ) );

    bool created = boost::filesystem::create_directory( testPath );

    EXPECT_TRUE( created );
    EXPECT_TRUE( boost::filesystem::is_directory( testPath ) );

    {
        scoped_file_cleanup fileCleanup;
        fileCleanup.add( testPath );
    }

    EXPECT_FALSE( boost::filesystem::exists( testPath ) );
}

TEST( ScopedFileCleanup, CleanFilledDirectory ) {
    using namespace frantic::files;

    boost::filesystem::path testPath = boost::filesystem::unique_path(
        boost::filesystem::temp_directory_path() / boost::filesystem::path( L"%%%%-%%%%-%%%%-%%%%" ) );

    EXPECT_FALSE( boost::filesystem::exists( testPath ) );

    bool created = boost::filesystem::create_directory( testPath );
    EXPECT_TRUE( created );
    EXPECT_TRUE( boost::filesystem::is_directory( testPath ) );

    for( size_t i = 0; i < 10; ++i ) {
        boost::filesystem::path tempFile =
            boost::filesystem::unique_path( testPath / boost::filesystem::path( "%%%%-%%%%-%%%%-%%%%.tmp" ) );
        std::ofstream os( testPath.c_str() );
        os << "Lorem ipsum";
        os.close();
    }

    {
        scoped_file_cleanup fileCleanup;
        fileCleanup.add( testPath );
    }

    EXPECT_FALSE( boost::filesystem::exists( testPath ) );
}
