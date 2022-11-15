// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#ifdef OIIO_LIB_AVAILABLE

#include <gtest/gtest.h>

#include <utilities/mesh_generators.hpp>

#include <frantic/geometry/trimesh3_interface.hpp>

#include <frantic/geometry/spatial_sampler.hpp>

#include <frantic/projection/texture_baking.hpp>

#include <frantic/graphics/color_rgba_f.hpp>
#include <frantic/graphics/vector3f.hpp>
#include <frantic/graphics2d/vector2f.hpp>

#include <OpenImageIO/imagebuf.h>

namespace {

class direct_image_projection : public frantic::geometry::spatial_sampler_channel {
  public:
    virtual frantic::channels::data_type_t data_type() const { return frantic::channels::data_type_float32; }

    virtual size_t arity() const { return 4; }

    virtual bool sample( const frantic::graphics::vector3fd& position, const frantic::graphics::vector3fd& /*normal*/,
                         char* outData ) const {
        using namespace frantic::graphics;
        using namespace frantic::graphics2d;

        color_rgba_f outPixel( float( position.x ), float( position.y ), float( position.z ), 1.0f );
        *reinterpret_cast<color_rgba_f*>( outData ) = outPixel;
        return true;
    }

    virtual bool has_density() const { return false; }
};

} // namespace

TEST( TextureBaking, SimpleBake ) {
    using namespace frantic::projection;
    using namespace frantic::graphics;
    using namespace frantic::graphics2d;
    using namespace frantic::geometry;
    using namespace OpenImageIO;

    direct_image_projection proj;

    trimesh3 mesh;
    make_quad_trimesh( mesh );
    mesh.add_vertex_channel<vector2f>( _T( "TextureCoord" ) );
    mesh.build_vertex_normals();

    trimesh3_vertex_channel_accessor<vector2f> texCoord(
        mesh.get_vertex_channel_accessor<frantic::graphics2d::vector2f>( _T( "TextureCoord" ) ) );

    for( size_t i = 0; i < 4; ++i ) {
        frantic::graphics2d::vector2f v = mesh.get_vertex( i ).project_xy();
        texCoord[i] = ( mesh.get_vertex( i ).project_xy() * 0.8f ) + vector2f( 0.1f, 0.1f );
    }

    mesh_interface_ptr iMesh( trimesh3_interface::create_instance( boost::move( mesh ) ).release() );

    ImageBuf image( ImageSpec( 256, 256, 4, TypeDesc::FLOAT ) );

    atlas_baking_params params;
    params.pixel_offset( vector2f( 0.0f, 0.0f ) ).overfill_distance_pixels( 3 );

    bake_atlas_texture( &proj, iMesh, image, params );

    image.write( "test.png", "png" );

    ImageSpec s;
    s.x = int( image.spec().width * 0.1f );
    s.y = int( image.spec().height * 0.1f );
    s.width = int( image.spec().width * 0.8f );
    s.height = int( image.spec().height * 0.8f );

    ROI roi( int( float( image.spec().width ) * 0.1f ), int( float( image.spec().width ) * 0.9f ),
             int( float( image.spec().height ) * 0.1f ), int( float( image.spec().height ) * 0.9f ) );

    typedef ImageBuf::Iterator<color_rgba_f, color_rgba_f> ImageIterator;

    for( ImageIterator it( image, roi ); !it.done(); ++it ) {
        color_rgba_f pixel = *( *it ).get();

        EXPECT_NEAR( float( it.x() - roi.xbegin ) / float( roi.xend - roi.xbegin ), pixel.get_r(), 0.01f );
        // Because the image space is flipped on the y-axis compared to texture space
        EXPECT_NEAR( float( it.y() - roi.ybegin ) / float( roi.yend - roi.ybegin ), ( 1.0f - pixel.get_g() ), 0.01f );
    }

    ROI pixelsOut( roi.xbegin - (int)params.m_overfillDistancePixels, roi.xend + (int)params.m_overfillDistancePixels,
                   roi.ybegin - (int)params.m_overfillDistancePixels, roi.yend + (int)params.m_overfillDistancePixels );

    for( ImageIterator it( image, pixelsOut ); !it.done(); ++it ) {
        EXPECT_LT( 0.0f, ( *it ).get()->get_a() );
    }
}

#endif
