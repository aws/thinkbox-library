// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/const_shared_polymesh3_builder.hpp>

#include <frantic/geometry/polymesh3.hpp>

namespace frantic {
namespace geometry {

const_shared_polymesh3_builder::const_shared_polymesh3_builder( polymesh3_channel_data& vertData,
                                                                polymesh3_channel_faces& geomPolyIndices,
                                                                polymesh3_channel_faces& geomPolyEndIndices ) {
    m_mesh = new polymesh3( vertData, geomPolyIndices, geomPolyEndIndices );
}

void const_shared_polymesh3_builder::add_vertex_channel( const frantic::tstring& channel,
                                                         polymesh3_channel_data& inVertexBuffer,
                                                         polymesh3_channel_faces* pInFaceBuffer ) {
    m_mesh->add_vertex_channel( channel, inVertexBuffer, pInFaceBuffer );
}

void const_shared_polymesh3_builder::add_face_channel( const frantic::tstring& channel,
                                                       polymesh3_channel_data& inFaceBuffer ) {
    m_mesh->add_face_channel( channel, inFaceBuffer );
}

const_polymesh3_ptr const_shared_polymesh3_builder::finalize() { return m_mesh; }

} // namespace geometry
} // namespace frantic
