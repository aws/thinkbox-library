// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <frantic/files/files.hpp>

using namespace std;
using namespace frantic;

TEST( FilesTest, PathManipulation ) {
#ifdef _WIN32
    std::string filename = "z:\\path\\to\\file.ext", numberedFilename = "z:\\path\\to\\file00343.tif";

    ASSERT_EQ( files::basename_from_path( filename ), "file" );
    ASSERT_EQ( files::directory_from_path( filename ), "z:\\path\\to" );
    ASSERT_EQ( files::ensure_trailing_pathseparator( filename ), "z:\\path\\to\\file.ext\\" );
    ASSERT_EQ( files::extension_from_path( filename ), ".ext" );
    ASSERT_EQ( files::filename_from_path( filename ), "file.ext" );
    ASSERT_EQ( files::forward_slashes( filename ), "z:/path/to/file.ext" );
    ASSERT_EQ( files::normalized_directory_name( "z:\\/path//to\\\\file.ext" ), filename );
    ASSERT_EQ( files::replace_directory( filename, "z:\\path" ), "z:\\path\\file.ext" );
    ASSERT_EQ( files::replace_extension( filename, ".tif" ), "z:\\path\\to\\file.tif" );
    ASSERT_EQ( files::replace_filename( filename, "other.tif" ), "z:\\path\\to\\other.tif" );
    ASSERT_EQ( files::replace_sequence_number( filename, 15 ), "z:\\path\\to\\file0015.ext" );
    ASSERT_EQ( files::replace_sequence_number( filename, 40015 ), "z:\\path\\to\\file40015.ext" );
    ASSERT_EQ( files::replace_sequence_number( filename, -15 ), "z:\\path\\to\\file-015.ext" );
    ASSERT_EQ( files::extract_sequence_number( numberedFilename ), 343 );
    ASSERT_EQ( files::increment_sequence_number( numberedFilename ), "z:\\path\\to\\file00344.tif" );

    std::string prefix, postfix;
    int framePadding, sequenceNumber;
    files::split_sequence_path( numberedFilename, prefix, framePadding, sequenceNumber, postfix );
    ASSERT_EQ( prefix, "z:\\path\\to\\file" );
    ASSERT_EQ( framePadding, 5 );
    ASSERT_EQ( sequenceNumber, 343 );
    ASSERT_EQ( postfix, ".tif" );

    ASSERT_EQ( files::directory_exists( filename ), false );
    ASSERT_EQ( files::file_exists( filename ), false );
    ASSERT_EQ( files::file_exists( "\\\\nonexistent_machine\\unc\\path.txt" ), false );
#else
    std::string filename = "/path/to/file.ext", numberedFilename = "/path/to/file00343.tif";

    ASSERT_EQ( files::basename_from_path( filename ), "file" );
    ASSERT_EQ( files::directory_from_path( filename ), "/path/to" );
    ASSERT_EQ( files::ensure_trailing_pathseparator( filename ), "/path/to/file.ext/" );
    ASSERT_EQ( files::extension_from_path( filename ), ".ext" );
    ASSERT_EQ( files::filename_from_path( filename ), "file.ext" );
    ASSERT_EQ( files::forward_slashes( filename ), "/path/to/file.ext" );
    ASSERT_EQ( files::normalized_directory_name( "//path//to//file.ext" ), filename );
    ASSERT_EQ( files::replace_directory( filename, "/path" ), "/path/file.ext" );
    ASSERT_EQ( files::replace_extension( filename, ".tif" ), "/path/to/file.tif" );
    ASSERT_EQ( files::replace_filename( filename, "other.tif" ), "/path/to/other.tif" );
    ASSERT_EQ( files::replace_sequence_number( filename, 15 ), "/path/to/file0015.ext" );
    ASSERT_EQ( files::replace_sequence_number( filename, 40015 ), "/path/to/file40015.ext" );
    ASSERT_EQ( files::replace_sequence_number( filename, -15 ), "/path/to/file-015.ext" );
    ASSERT_EQ( files::extract_sequence_number( numberedFilename ), 343 );
    ASSERT_EQ( files::increment_sequence_number( numberedFilename ), "/path/to/file00344.tif" );

    std::string prefix, postfix;
    int framePadding, sequenceNumber;
    files::split_sequence_path( numberedFilename, prefix, framePadding, sequenceNumber, postfix );
    ASSERT_EQ( prefix, "/path/to/file" );
    ASSERT_EQ( framePadding, 5 );
    ASSERT_EQ( sequenceNumber, 343 );
    ASSERT_EQ( postfix, ".tif" );

    ASSERT_EQ( files::directory_exists( filename ), false );
    ASSERT_EQ( files::file_exists( filename ), false );
#endif
}

TEST( FilesTest, MakeRelativePath ) {
#ifdef _WIN32
    string root = "C:/root/Path", path = "C:\\Root//path\\some\\File.txt";
    EXPECT_TRUE( files::make_relative_path( root, path ) );
    EXPECT_EQ( "some\\File.txt", path );

    path = "C:/root\\Other/File.txt";
    EXPECT_TRUE( files::make_relative_path( root, path ) );
    EXPECT_EQ( "..\\Other\\File.txt", path );

    path = "D:/root/Path/File.txt";
    EXPECT_FALSE( files::make_relative_path( root, path ) );
    EXPECT_EQ( "D:/root/Path/File.txt", path );
#else
    string root = "/root/Path", path = "/root//Path/some/File.txt";
    EXPECT_TRUE( files::make_relative_path( root, path ) );
    EXPECT_EQ( "some/File.txt", path );

    path = "/toot/Path/File.txt";
    EXPECT_TRUE( files::make_relative_path( root, path ) );
    EXPECT_EQ( "../../toot/Path/File.txt", path );
#endif
}

TEST( FilesTest, ExtensionCompare ) {
    EXPECT_TRUE( files::extension_iequals( "files.ext", ".ext" ) );
    EXPECT_TRUE( files::extension_iequals( "file.ExT", ".ext" ) );
    EXPECT_TRUE( files::extension_iequals( "file.ExT", ".eXt" ) );
    EXPECT_FALSE( files::extension_iequals( "file.exs", ".ext" ) );
}
