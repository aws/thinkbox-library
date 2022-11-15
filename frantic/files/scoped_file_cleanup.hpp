// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/strings/tstring.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/move/move.hpp>

#include <vector>

namespace frantic {
namespace files {

/**
 * A simple RAII class to ensure temporary files are deleted when a function goes
 * out of scope. A replacemen
 *
 * To consider:
 * - make this movable-but-non-copyable
 * - make this mergable (allow a number of tasks to create temp files, and then delete them all when the operation is
 * complete)
 */
class scoped_file_cleanup {
  private:
    // for now at least, non-copyable
    scoped_file_cleanup( const scoped_file_cleanup& );
    scoped_file_cleanup operator=( const scoped_file_cleanup& );

  public:
    scoped_file_cleanup();
    ~scoped_file_cleanup();

    void reset();

    void add( const frantic::tstring& filePath );
    void add( const boost::filesystem::path& filePath );

    void cleanup_files();

  private:
    std::vector<boost::filesystem::path> m_deletionPaths;
};

} // namespace files
} // namespace frantic
