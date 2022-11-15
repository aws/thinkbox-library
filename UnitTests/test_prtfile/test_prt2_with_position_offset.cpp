// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "tbb/task_scheduler_init.h"

#include <frantic/particles/prt_metadata.hpp>
#include <frantic/prtfile/prt2_reader.hpp>
#include <frantic/prtfile/prt2_writer.hpp>

#pragma warning( push )
#pragma warning( disable : 4100 4512 )
#include <tbb/pipeline.h>
#pragma warning( pop )

using namespace std;
using frantic::channels::channel_cvt_accessor;
using frantic::channels::channel_map;
using frantic::channels::property_map;
using frantic::graphics::boundbox3f;
using frantic::graphics::color3h;
using frantic::graphics::vector3f;
using frantic::particles::particle_array;
using frantic::prtfile::prt2_reader;
using frantic::prtfile::prt2_writer;

// A handcrafted PRT2 file with particle layout {Position: 3 * float32, Color: 3 * float16, Density: float32},
// 8 particles in 3 chunks, and with some global and channel metadata.
// This example is using a 'PrtO' filechunk for the particles instead of 'Part'.
//
// clang-format off
static unsigned const char prtdata[] = {
  // Magic value
  0xC0, 'P', 'R', 'T', '2', 0x0D, 0x0A, 0x1A,
  // PRT version number
  0x03, 0x00, 0x00, 0x00,

  // >> FileChunk 'Chan'
  'C', 'h', 'a', 'n',
  //    chunkSize (uint64)
  0x3B, 0, 0, 0, 0, 0, 0, 0,
  //    chunkData
  //      channelCount (varint)
  3,
  //        channelName (varstring)
  8,  'P', 'o', 's', 'i', 't', 'i', 'o', 'n',
  //        dataType (varstring)
  11, '3', ' ', '*', ' ', 'f', 'l', 'o', 'a', 't', '3',
      '2',
  //        sizeBytes (varint)
  12,
  //        channelName (varstring)
  5,  'C', 'o', 'l', 'o', 'r',
  //        dataType (varstring)
  11, '3', ' ', '*', ' ', 'f', 'l', 'o', 'a', 't', '1',
      '6',
  //        sizeBytes (varint)
  6,
  //        channelName (varstring)
  7,  'D', 'e', 'n', 's', 'i', 't', 'y',
  //        dataType (varstring)
  7,  'f', 'l', 'o', 'a', 't', '3', '2',
  //        sizeBytes (varint)
  4,

  // >> FileChunk 'PrtO'
  'P', 'r', 't', 'O',
  //    chunkSize (uint64)
  0x0a, 0x01, 0, 0, 0, 0, 0, 0,
  //    particleStreamName (varstring), "" = default stream
  0,
  //    compressionScheme (varstring), "uncompressed" = no compression
  12, 'u', 'n', 'c', 'o', 'm', 'p', 'r', 'e', 's', 's',
      'e', 'd',
  //    particleCount (uint64), 8 particles
  0x08, 0, 0, 0, 0, 0, 0, 0,
  //    particleChunkCount (uint64), 3 chunks
  0x03, 0, 0, 0, 0, 0, 0, 0,
  //      chunk 0 size in bytes (uint32)
  56, 0, 0, 0,
  //      chunk 0 particle count (uint32)
  2, 0, 0, 0,
  //      chunk 0 PositionOffset (1.5, -1, 2)
  0x00, 0x00, 0xc0, 0x3f,  0x00, 0x00, 0x80, 0xBF,  0x00, 0x00, 0x00, 0x40,
  //      0 Position
  0x00, 0x00, 0x80, 0xBF,  0x00, 0x00, 0x80, 0xBF,  0x00, 0x00, 0x80, 0xBF,
  //      0 Color
  0x00, 0x38,  0x00, 0x3c,  0x00, 0x38,
  //      0 Density
  0x00, 0x00, 0x80, 0x3F,
  //      1 Position
  0x00, 0x00, 0x80, 0xBF,  0x00, 0x00, 0x80, 0xBF,  0x00, 0x00, 0x80, 0x3F,
  //      1 Color
  0x00, 0x38,  0x00, 0x3c,  0x00, 0x38,
  //      1 Density
  0x00, 0x00, 0x80, 0x3F,
  //      chunk 1 size in bytes (uint32)
  78, 0, 0, 0,
  //      chunk 1 particle count (uint32)
  3, 0, 0, 0,
  //      chunk 1 PositionOffset (0, 0, 0)
  0x00, 0x00, 0x00, 0x00,  0x00, 0x00, 0x00, 0x00,  0x00, 0x00, 0x00, 0x00,
  //      2 Position
  0x00, 0x00, 0x80, 0xBF,  0x00, 0x00, 0x80, 0x3F,  0x00, 0x00, 0x80, 0xBF,
  //      2 Color
  0x00, 0x38,  0x00, 0x3c,  0x00, 0x38,
  //      2 Density
  0x00, 0x00, 0x80, 0x3F,
  //      3 Position
  0x00, 0x00, 0x80, 0xBF,  0x00, 0x00, 0x80, 0x3F,  0x00, 0x00, 0x80, 0x3F,
  //      3 Color
  0x00, 0x38,  0x00, 0x3c,  0x00, 0x38,
  //      3 Density
  0x00, 0x00, 0x80, 0x3F,
  //      4 Position
  0x00, 0x00, 0x80, 0x3F,  0x00, 0x00, 0x80, 0xBF,  0x00, 0x00, 0x80, 0xBF,
  //      4 Color
  0x00, 0x38,  0x00, 0x3c,  0x00, 0x38,
  //      4 Density
  0x00, 0x00, 0x80, 0x3F,
  //      chunk 2 size in bytes (uint32)
  78, 0, 0, 0,
  //      chunk 2 particle count (uint32)
  3, 0, 0, 0,
  //      chunk 2 PositionOffset (1, 2, 3)
  0x00, 0x00, 0x80, 0x3f,  0x00, 0x00, 0x00, 0x40,  0x00, 0x00, 0x40, 0x40,
  //      5 Position
  0x00, 0x00, 0x80, 0x3F,  0x00, 0x00, 0x80, 0xBF,  0x00, 0x00, 0x80, 0x3F,
  //      5 Color
  0x00, 0x38,  0x00, 0x3c,  0x00, 0x38,
  //      5 Density
  0x00, 0x00, 0x80, 0x3F,
  //      6 Position
  0x00, 0x00, 0x80, 0x3F,  0x00, 0x00, 0x80, 0x3F,  0x00, 0x00, 0x80, 0xBF,
  //      6 Color
  0x00, 0x38,  0x00, 0x3c,  0x00, 0x38,
  //      6 Density
  0x00, 0x00, 0x80, 0x3F,
  //      7 Position
  0x00, 0x00, 0x80, 0x3F,  0x00, 0x00, 0x80, 0x3F,  0x00, 0x00, 0x80, 0x3F,
  //      7 Color
  0x00, 0x38,  0x00, 0x3c,  0x00, 0x38,
  //      7 Density
  0x00, 0x00, 0x80, 0x3F,

  // >> FileChunk 'PIdx'
  'P', 'I', 'd', 'x',
  //    chunkSize (uint64)
  0x0f, 0, 0, 0, 0, 0, 0, 0,
  //    particleStreamName (varstring), "" = default stream
  0,
  //    particleChunkCount (uint64), 2 chunks
  0x03, 0, 0, 0, 0, 0, 0, 0,
  //    chunk 0 size
  64,
  //    chunk 0 particleCount
  2,
  //    chunk 1 size
  86,
  //    chunk 1 particleCount
  3,
  //    chunk 2 size
  86,
  //    chunk 2 particleCount
  3,
};
// clang-format on

TEST( PRT2WithPositionOffsets, ReadHandcraftedFile ) {
    prt2_reader prt2;
    prt2.open( new istringstream( string( prtdata, prtdata + sizeof( prtdata ) ) ), _T("handcrafted.prt") );

    // Particle/chunk counts
    EXPECT_EQ( 8, prt2.get_particle_count() );
    EXPECT_EQ( 3, prt2.get_particle_chunk_count() );

    // Verify the particle chunk index
    pair<boost::int64_t, boost::int64_t> ex;
    ex = prt2.get_particle_chunk_extents( 0 );
    EXPECT_EQ( 0, ex.first );
    EXPECT_EQ( 2, ex.second );
    ex = prt2.get_particle_chunk_extents( 1 );
    EXPECT_EQ( 2, ex.first );
    EXPECT_EQ( 5, ex.second );
    ex = prt2.get_particle_chunk_extents( 2 );
    EXPECT_EQ( 5, ex.first );
    EXPECT_EQ( 8, ex.second );

    particle_array parr( prt2.get_channel_map() );
    channel_cvt_accessor<vector3f> posAccessor = prt2.get_channel_map().get_cvt_accessor<vector3f>( _T("Position") );

    // Check the contents of the particles themselves, allowing the prt2_reader to apply the offset
    prt2.read_particle_chunk( 0, parr );
    EXPECT_EQ( 2u, parr.particle_count() );
    EXPECT_EQ( vector3f( -1, -1, -1 ) + vector3f( 1.5f, -1.f, 2.f ), posAccessor( parr[0] ) );
    EXPECT_EQ( vector3f( -1, -1, 1 ) + vector3f( 1.5f, -1.f, 2.f ), posAccessor( parr[1] ) );
    prt2.read_particle_chunk( 1, parr );
    EXPECT_EQ( 3u, parr.particle_count() );
    EXPECT_EQ( vector3f( -1, 1, -1 ), posAccessor( parr[0] ) );
    EXPECT_EQ( vector3f( -1, 1, 1 ), posAccessor( parr[1] ) );
    EXPECT_EQ( vector3f( 1, -1, -1 ), posAccessor( parr[2] ) );
    prt2.read_particle_chunk( 2, parr );
    EXPECT_EQ( 3u, parr.particle_count() );
    EXPECT_EQ( vector3f( 1, -1, 1 ) + vector3f( 1.f, 2.f, 3.f ), posAccessor( parr[0] ) );
    EXPECT_EQ( vector3f( 1, 1, -1 ) + vector3f( 1.f, 2.f, 3.f ), posAccessor( parr[1] ) );
    EXPECT_EQ( vector3f( 1, 1, 1 ) + vector3f( 1.f, 2.f, 3.f ), posAccessor( parr[2] ) );

    // Check the contents of the particles themselves, also retrieving the offsets separately
    vector3f positionOffset;
    prt2.read_particle_chunk( 0, parr, &positionOffset );
    EXPECT_EQ( 2u, parr.particle_count() );
    EXPECT_EQ( vector3f( 1.5f, -1.f, 2.f ), positionOffset );
    EXPECT_EQ( vector3f( -1, -1, -1 ), posAccessor( parr[0] ) );
    EXPECT_EQ( vector3f( -1, -1, 1 ), posAccessor( parr[1] ) );
    prt2.read_particle_chunk( 1, parr, &positionOffset );
    EXPECT_EQ( 3u, parr.particle_count() );
    EXPECT_EQ( vector3f( 0.f ), positionOffset );
    EXPECT_EQ( vector3f( -1, 1, -1 ), posAccessor( parr[0] ) );
    EXPECT_EQ( vector3f( -1, 1, 1 ), posAccessor( parr[1] ) );
    EXPECT_EQ( vector3f( 1, -1, -1 ), posAccessor( parr[2] ) );
    prt2.read_particle_chunk( 2, parr, &positionOffset );
    EXPECT_EQ( vector3f( 1.f, 2.f, 3.f ), positionOffset );
    EXPECT_EQ( 3u, parr.particle_count() );
    EXPECT_EQ( vector3f( 1, -1, 1 ), posAccessor( parr[0] ) );
    EXPECT_EQ( vector3f( 1, 1, -1 ), posAccessor( parr[1] ) );
    EXPECT_EQ( vector3f( 1, 1, 1 ), posAccessor( parr[2] ) );
}

TEST( PRT2WithPositionOffsets, WriteHandcraftedFile ) {
    tbb::task_scheduler_init taskScheduleInit;
    channel_map pcm;
    pcm.end_channel_definition();
    pcm.append_channel<vector3f>( _T("Position") );
    pcm.append_channel<color3h>( _T("Color") );
    pcm.append_channel<float>( _T("Density") );

    // Write all the data that's in the handcrafted PRT2 file, with chunks in the same order, etc
    prt2_writer prt2;
    ostringstream* os = new ostringstream();
    prt2.open( os, _T("handcrafted.prt"), pcm );

    particle_array parr( prt2.get_channel_map() );
    parr.resize( 8 );
    channel_cvt_accessor<vector3f> pos = prt2.get_channel_map().get_cvt_accessor<vector3f>( _T("Position") );
    channel_cvt_accessor<color3h> col = prt2.get_channel_map().get_cvt_accessor<color3h>( _T("Color") );
    channel_cvt_accessor<float> dens = prt2.get_channel_map().get_cvt_accessor<float>( _T("Density") );
    for( int i = 0; i < 8; ++i ) {
        pos.set( parr.at( i ), vector3f( ( i & 4 ) ? 1.f : -1.f, ( i & 2 ) ? 1.f : -1.f, ( i & 1 ) ? 1.f : -1.f ) );
        col.set( parr.at( i ), color3h( 0.5f, 1.f, 0.5f ) );
        dens.set( parr.at( i ), 1.f );
    }

    // A chunk generator that takes the handcrafted data and allows it to be accessed in a pull format.
    class chunk_gen : public tbb::filter {
        particle_array& m_particles;
        int m_i;

      public:
        chunk_gen( particle_array& particles )
            : tbb::filter( true )
            , m_particles( particles )
            , m_i( -1 ) {}

        void* operator()( void* ) {
            switch( ++m_i ) {
            case 0: {
                prt2_writer::particle_chunk* chunk = new prt2_writer::particle_chunk();
                chunk->positionOffset.set( 1.5f, -1.f, 2.f );
                chunk->usePositionOffset = true;
                chunk->particleCount = 2;
                chunk->uncompressed.resize( m_particles.get_channel_map().structure_size() * 2 );
                memcpy( &chunk->uncompressed[0], m_particles.at( 0 ), chunk->uncompressed.size() );
                return chunk;
            }
            case 1: {
                prt2_writer::particle_chunk* chunk = new prt2_writer::particle_chunk();
                chunk->positionOffset.set( 0.f );
                chunk->usePositionOffset = true;
                chunk->particleCount = 3;
                chunk->uncompressed.resize( m_particles.get_channel_map().structure_size() * 3 );
                memcpy( &chunk->uncompressed[0], m_particles.at( 2 ), chunk->uncompressed.size() );
                return chunk;
            }
            case 2: {
                prt2_writer::particle_chunk* chunk = new prt2_writer::particle_chunk();
                chunk->positionOffset.set( 1, 2, 3 );
                chunk->usePositionOffset = true;
                chunk->particleCount = 3;
                chunk->uncompressed.resize( m_particles.get_channel_map().structure_size() * 3 );
                memcpy( &chunk->uncompressed[0], m_particles.at( 5 ), chunk->uncompressed.size() );
                return chunk;
            }
            case 3: {
                return &prt2_writer::TERMINATION_CHUNK;
            }
            default: {
                return NULL;
            }
            }
        }
    };

    boost::shared_ptr<tbb::filter> chunkGen( new chunk_gen( parr ) );
    frantic::logging::null_progress_logger nullProgress;
    prt2.write_particle_chunks( chunkGen, 8, nullProgress, true, frantic::prtfile::prt2_compression_uncompressed );

    // Compare what we've written in the stream with the handcrafted data
    os->flush();
    string s = os->str();
    prt2.close();

    /*
    // If debugging is necessary, use a hexdump tool like HxD to compare the two files:
    ofstream fout( "x.prt", ios::out | ios::binary );
    fout.write( s.data(), s.size() );
    fout.close();
    fout.open( "y.prt", ios::out | ios::binary );
    fout.write( (const char *)prtdata, sizeof( prtdata ) );
    //*/
    ASSERT_EQ( sizeof( prtdata ), s.size() );
    EXPECT_TRUE( memcmp( prtdata, s.data(), s.size() ) == 0 );
}
