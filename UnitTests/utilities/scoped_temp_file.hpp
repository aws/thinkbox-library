// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/filesystem.hpp>

#include <frantic/files/files.hpp>

class scoped_temp_file : boost::noncopyable {
  public:
    scoped_temp_file( boost::filesystem::path extension )
        : m_path( boost::filesystem::unique_path(
              ( boost::filesystem::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%." ).replace_extension( extension ) ) ) {
    }

    ~scoped_temp_file() {
        boost::system::error_code ec;
        boost::filesystem::remove( m_path, ec );
    }

    frantic::tstring get_path() const { return frantic::files::to_tstring( m_path ); }

  private:
    boost::filesystem::path m_path;
};