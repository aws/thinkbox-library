// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/files/scoped_file_cleanup.hpp>

#include <boost/filesystem/operations.hpp>

namespace frantic {
namespace files {

scoped_file_cleanup::scoped_file_cleanup() {}

scoped_file_cleanup::~scoped_file_cleanup() { cleanup_files(); }

void scoped_file_cleanup::reset() { m_deletionPaths.clear(); }

void scoped_file_cleanup::add( const frantic::tstring& filePath ) { add( boost::filesystem::path( filePath ) ); }

void scoped_file_cleanup::add( const boost::filesystem::path& filePath ) { m_deletionPaths.push_back( filePath ); }

void scoped_file_cleanup::cleanup_files() {
    for( size_t i = 0; i < m_deletionPaths.size(); ++i ) {
        if( boost::filesystem::exists( m_deletionPaths[i] ) ) {
            boost::filesystem::remove_all( m_deletionPaths[i] );
        }
    }

    m_deletionPaths.clear();
}

} // namespace files
} // namespace frantic
