// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include <stdafx.h>
// clang-format on

#include <frantic/channels/channel_map.hpp>
#include <frantic/geometry/triangle_utils.hpp>
#include <frantic/particles/streams/surface_particle_istream.hpp>

using namespace frantic::particles::streams;

namespace {
/**
 * Wrapper used for creation of surface_particle_istream.
 */
class trimesh3_wrapper {
    boost::shared_ptr<frantic::geometry::trimesh3> m_mesh;
    frantic::channels::channel_accessor<frantic::graphics::vector3f> m_posAccessor;
    frantic::channels::channel_cvt_accessor<frantic::graphics::vector3f> m_normalAccessor;

  public:
    trimesh3_wrapper( boost::shared_ptr<frantic::geometry::trimesh3> mesh ) { m_mesh = mesh; }

    ~trimesh3_wrapper() {}

    void get_native_map( frantic::channels::channel_map& outNativeMap ) {
        outNativeMap.define_channel<frantic::graphics::vector3f>( _T("Position") );
        outNativeMap.define_channel<frantic::graphics::vector3f>( _T("Normal") );
    }

    void set_channel_map( const frantic::channels::channel_map& seedMap ) {
        m_normalAccessor.reset();
        m_posAccessor.reset();

        if( seedMap.has_channel( _T( "Position" ) ) )
            m_posAccessor = seedMap.get_accessor<frantic::graphics::vector3f>( _T( "Position" ) );

        if( seedMap.has_channel( _T( "Normal" ) ) )
            m_normalAccessor = seedMap.get_cvt_accessor<frantic::graphics::vector3f>( _T( "Normal" ) );
    }

    std::size_t surface_count() { return m_mesh->face_count(); }

    std::size_t element_count( std::size_t surfaceIndex ) {
        if( 0 <= surfaceIndex && surfaceIndex < m_mesh->face_count() )
            return 1; // Triangles is only split into 1 triangle
        return 0;
    }

    float element_area( std::size_t surfaceIndex, std::size_t /*elementIndex*/ ) {
        float area = -1.0f;
        frantic::geometry::vector3 vertices;
        frantic::geometry::vector3f vertex1;
        frantic::geometry::vector3f vertex2;
        frantic::geometry::vector3f vertex3;

        if( 0 <= surfaceIndex && surfaceIndex < m_mesh->face_count() ) {
            vertices = m_mesh->get_face( surfaceIndex );

            vertex1 = m_mesh->get_vertex( vertices.x );
            vertex2 = m_mesh->get_vertex( vertices.y );
            vertex3 = m_mesh->get_vertex( vertices.z );

            area = 0.5f * frantic::geometry::vector3f::cross( vertex1 - vertex2, vertex3 - vertex2 ).get_magnitude();
        }

        return area;
    }

    template <class RandomGen>
    void seed_particle( char* pOutParticle, std::size_t surfaceIndex, std::size_t /*elementIndex*/,
                        RandomGen& randomnessGenerator ) {
        float baryCentricCoord[3];
        frantic::graphics::vector3f triVerts[3];
        frantic::geometry::vector3 vertices;
        frantic::geometry::vector3f particle;

        frantic::geometry::random_barycentric_coordinate( baryCentricCoord, randomnessGenerator );

        vertices = m_mesh->get_face( surfaceIndex );

        triVerts[0] = m_mesh->get_vertex( vertices.x );
        triVerts[1] = m_mesh->get_vertex( vertices.y );
        triVerts[2] = m_mesh->get_vertex( vertices.z );

        particle = ( baryCentricCoord[0] * triVerts[0] + baryCentricCoord[1] * triVerts[1] +
                     baryCentricCoord[2] * triVerts[2] );

        m_posAccessor.get( pOutParticle ) = particle;

        if( m_normalAccessor.is_valid() )
            m_normalAccessor.set( pOutParticle,
                                  frantic::graphics::triangle_normal( triVerts[0], triVerts[1], triVerts[2] ) );
    }
};

} // anonymous namespace

particle_istream_ptr frantic::particles::streams::create_surface_particle_istream_using_count(
    const frantic::geometry::trimesh3& mesh, const boost::uint64_t particleCount, const boost::uint32_t randomSeed ) {

    boost::shared_ptr<frantic::geometry::trimesh3> meshPtr( new frantic::geometry::trimesh3( mesh ) );
    trimesh3_wrapper trimeshWrapper( meshPtr );
    boost::shared_ptr<surface_particle_istream<trimesh3_wrapper>> surfaceStream(
        new surface_particle_istream<trimesh3_wrapper>( trimeshWrapper ) );

    surfaceStream->set_random_seed( randomSeed );
    surfaceStream->set_particle_count( particleCount );

    return surfaceStream;
}

particle_istream_ptr frantic::particles::streams::create_surface_particle_istream_using_spacing(
    const frantic::geometry::trimesh3& mesh, const float particleSpacing, const boost::uint32_t randomSeed ) {

    boost::shared_ptr<frantic::geometry::trimesh3> meshPtr( new frantic::geometry::trimesh3( mesh ) );
    trimesh3_wrapper trimeshWrapper( meshPtr );
    boost::shared_ptr<surface_particle_istream<trimesh3_wrapper>> surfaceStream(
        new surface_particle_istream<trimesh3_wrapper>( trimeshWrapper ) );

    surfaceStream->set_random_seed( randomSeed );
    surfaceStream->set_particle_spacing( particleSpacing );

    return surfaceStream;
}
