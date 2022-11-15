// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/geometry/polymesh3_builder.hpp>

namespace frantic {
namespace geometry {

void polymesh3_builder::removeNonPlanarPolygons() {
    throw std::runtime_error( "polymesh3_builder::removeNonPlanarPolygons() Not implemented yet" );
}

void polymesh3_builder::removeCoincidentVertices() {
    throw std::runtime_error( "polymesh3_builder::removeCoincidentVertices() Not implemented yet" );
}

void polymesh3_builder::removePolygonLoops() {
    throw std::runtime_error( "polymesh3_builder::removePolygonLoops() Not implemented yet" );
}

polymesh3_builder::polymesh3_builder() {}

polymesh3_ptr polymesh3_builder::finalize( bool removeNonPlanar, bool removeCoincident, bool removePolyRepeats ) {
    if( removeNonPlanar )
        removeNonPlanarPolygons();
    if( removeCoincident )
        removeCoincidentVertices();
    if( removePolyRepeats )
        removePolygonLoops();

    return new polymesh3( m_vertexBuffer, m_faceBuffer, m_faceEndOffsets );
}

void polymesh3_builder::add_vertex( float x, float y, float z ) {
    float* pVert = (float*)m_vertexBuffer.add_element( sizeof( float ) * 3 );
    pVert[0] = x;
    pVert[1] = y;
    pVert[2] = z;
}

void polymesh3_builder::add_vertex( const frantic::graphics::vector3f& p ) {
    float* pVert = (float*)m_vertexBuffer.add_element( sizeof( float ) * 3 );
    pVert[0] = p.x;
    pVert[1] = p.y;
    pVert[2] = p.z;
}

void polymesh3_builder::add_polygon( const std::vector<frantic::graphics::vector3f>& pts ) {
    add_polygon( &pts[0], pts.size() );
}

void polymesh3_builder::add_polygon( const frantic::graphics::vector3f* pPts, std::size_t count ) {
    std::size_t offset = m_vertexBuffer.size() / ( sizeof( float[3] ) );

    for( std::size_t i = 0; i < count; ++i ) {
        add_vertex( pPts[i] );
        m_faceBuffer.push_back( (int)( offset + i ) );
    }
    m_faceEndOffsets.push_back( (int)m_faceBuffer.size() );
}

void polymesh3_builder::add_polygon( const float ( *pPts )[3], std::size_t count ) {
    std::size_t offset = m_vertexBuffer.size() / ( sizeof( float[3] ) );

    for( std::size_t i = 0; i < count; ++i ) {
        add_vertex( pPts[i][0], pPts[i][1], pPts[i][2] );
        m_faceBuffer.push_back( (int)( offset + i ) );
    }
    m_faceEndOffsets.push_back( (int)m_faceBuffer.size() );
}

void polymesh3_builder::add_polygon( const std::vector<int>& ptIndices ) {
    add_polygon( &ptIndices[0], ptIndices.size() );
}

void polymesh3_builder::add_polygon( const int newIndices[], std::size_t count ) {
    for( std::size_t i = 0; i < count; ++i )
        m_faceBuffer.push_back( newIndices[i] );
    m_faceEndOffsets.push_back( (int)m_faceBuffer.size() );
}

} // namespace geometry
} // namespace frantic
