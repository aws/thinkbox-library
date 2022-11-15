// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/geometry/polymesh3.hpp>

namespace frantic {
namespace geometry {

// forward declaration
class ply_reader_header;

class ply_reader : boost::noncopyable {
  public:
    typedef std::pair<frantic::channels::data_type_t, std::size_t> channel_type_t;

    ply_reader( const frantic::tstring& filename );

    const frantic::tstring& get_filename() const;

    std::size_t get_vertex_count() const;
    std::size_t get_face_count() const;

    void get_vertex_channel_names( std::vector<frantic::tstring>& outNames ) const;

    channel_type_t get_vertex_channel_type( const frantic::tstring& name ) const;

    frantic::geometry::polymesh3_ptr read_polymesh3() const;

  private:
    ply_reader_header& get_header() const;

    frantic::tstring m_filename;

    boost::shared_ptr<ply_reader_header> m_header;
};

} // namespace geometry
} // namespace frantic
