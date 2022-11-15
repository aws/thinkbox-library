// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <fstream>

#include <tbb/task_scheduler_init.h>

#include <frantic/particles/prt_metadata.hpp>
#include <frantic/prtfile/prt2_reader.hpp>
#include <frantic/prtfile/prt2_writer.hpp>

#include <frantic/files/paths.hpp>
#include <frantic/files/scoped_file_cleanup.hpp>

#include <vector>

#pragma warning( push )
#pragma warning( disable : 4100 4512 )
#include <tbb/pipeline.h>
#pragma warning( pop )

using namespace std;
using namespace boost;
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
// 8 particles in 2 chunks, and with some global and channel metadata.
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

  // >> FileChunk 'Meta'  LengthUnitInMicrometers
  'M', 'e', 't', 'a',
  //    chunkSize (uint64)
  0x28, 0, 0, 0, 0, 0, 0, 0,
  //    name (varstring)
  23, 'L', 'e', 'n', 'g', 't', 'h', 'U', 'n', 'i', 't',
      'I', 'n', 'M', 'i', 'c', 'r', 'o', 'm', 'e', 't',
      'e', 'r', 's',
  //    type (varstring)
  7,  'f', 'l', 'o', 'a', 't', '6', '4',
  //    value (according to 'type'), float64 25400.0
  0x00, 0x00, 0x00, 0x00, 0x00, 0xce, 0xd8, 0x40,

  // >> FileChunk 'Meta'  CoordSys
  'M', 'e', 't', 'a',
  //    chunkSize (uint64)
  0x13, 0, 0, 0, 0, 0, 0, 0,
  //    name (varstring)
  8,  'C', 'o', 'o', 'r', 'd', 'S', 'y', 's',
  //    type (varstring)
  5,  'i', 'n', 't', '3', '2',
  //    value (according to 'type'), int32 2 right-handed Z-up
  0x02, 0x00, 0x00, 0x00,

  // >> FileChunk 'Meta'  FrameRate
  'M', 'e', 't', 'a',
  //    chunkSize (uint64)
  0x1d, 0, 0, 0, 0, 0, 0, 0,
  //    name (varstring)
  9,  'F', 'r', 'a', 'm', 'e', 'R', 'a', 't', 'e',
  //    type (varstring)
  10, '2', ' ', '*', ' ', 'u', 'i', 'n', 't', '3', '2',
  //    value (according to 'type'), uint32 pair (30000, 1001) = 29.97 NTSC framerate
  0x30, 0x75, 0x00, 0x00,
  0xE9, 0x03, 0x00, 0x00,

  // >> FileChunk 'Meta'  Position.Interpretation
  'M', 'e', 't', 'a',
  //    chunkSize (uint64)
  0x25, 0, 0, 0, 0, 0, 0, 0,
  //    name (varstring)
  23, 'P', 'o', 's', 'i', 't', 'i', 'o', 'n', '.', 'I',
      'n', 't', 'e', 'r', 'p', 'r', 'e', 't', 'a', 't',
      'i', 'o', 'n',
  //    type (varstring)
  6,  's', 't', 'r', 'i', 'n', 'g',
  //    value (according to 'type'), string
  5,  'P', 'o', 'i', 'n', 't',

  // >> FileChunk 'Part'
  'P', 'a', 'r', 't',
  //    chunkSize (uint64)
  0xde, 0, 0, 0, 0, 0, 0, 0,
  //    particleStreamName (varstring), "" = default stream
  0,
  //    compressionScheme (varstring), "uncompressed" = no compression
  12, 'u', 'n', 'c', 'o', 'm', 'p', 'r', 'e', 's', 's',
      'e', 'd',
  //    particleCount (uint64), 8 particles
  0x08, 0, 0, 0, 0, 0, 0, 0,
  //    particleChunkCount (uint64), 2 chunks
  0x02, 0, 0, 0, 0, 0, 0, 0,
  //      chunk 0 size in bytes (uint32)
  44, 0, 0, 0,
  //      chunk 0 particle count (uint32)
  2, 0, 0, 0,
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
  132, 0, 0, 0,
  //      chunk 1 particle count (uint32)
  6, 0, 0, 0,
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
  0x0e, 0, 0, 0, 0, 0, 0, 0,
  //    particleStreamName (varstring), "" = default stream
  0,
  //    particleChunkCount (uint64), 2 chunks
  0x02, 0, 0, 0, 0, 0, 0, 0,
  //    chunk 0 size
  52,
  //    chunk 0 particleCount
  2,
  //    chunk 1 size
  0x8c, 0x01,
  //    chunk 1 particleCount
  6,

  // >> FileChunk 'Meta'  Position.Extents
  'M', 'e', 't', 'a',
  //    chunkSize (uint64)
  0x35, 0, 0, 0, 0, 0, 0, 0,
  //    name (varstring)
  16, 'P', 'o', 's', 'i', 't', 'i', 'o', 'n', '.', 'E',
      'x', 't', 'e', 'n', 't', 's',
  //    type (varstring)
  11, '6', ' ', '*', ' ', 'f', 'l', 'o', 'a', 't', '3',
      '2',
  //    value (according to 'type'), xmin, ymin, zmin, xmax, ymax, zmax
  0x00, 0x00, 0x80, 0xBF, // -1.f
  0x00, 0x00, 0x80, 0xBF, // -1.f
  0x00, 0x00, 0x80, 0xBF, // -1.f
  0x00, 0x00, 0x80, 0x3F, // 1.f
  0x00, 0x00, 0x80, 0x3F, // 1.f
  0x00, 0x00, 0x80, 0x3F, // 1.f
};
// clang-format on

TEST( PRT2, ReadHandcraftedFile ) {
    prt2_reader prt2;
    prt2.open( new istringstream( string( prtdata, prtdata + sizeof( prtdata ) ) ), _T("handcrafted.prt") );

    // Particle/chunk counts
    EXPECT_EQ( 8, prt2.get_particle_count() );
    EXPECT_EQ( 2, prt2.get_particle_chunk_count() );

    // The channel map
    const channel_map& cm = prt2.get_channel_map();
    EXPECT_TRUE( cm.channel_definition_complete() );
    EXPECT_EQ( 3u, cm.channel_count() );
    EXPECT_EQ( 22u, cm.structure_size() );
    EXPECT_TRUE( cm.has_channel( _T("Position") ) );
    EXPECT_EQ( 0u, cm.channel_offset( _T("Position") ) );
    EXPECT_TRUE( cm.has_channel( _T("Color") ) );
    EXPECT_EQ( 12u, cm.channel_offset( _T("Color") ) );
    EXPECT_TRUE( cm.has_channel( _T("Density") ) );
    EXPECT_EQ( 18u, cm.channel_offset( _T("Density") ) );
    frantic::tstring name;
    frantic::channels::data_type_t dtype;
    size_t arity;
    cm.get_channel_definition( 0, name, dtype, arity );
    EXPECT_EQ( _T("Position"), name );
    EXPECT_EQ( frantic::channels::data_type_float32, dtype );
    EXPECT_EQ( 3u, arity );
    cm.get_channel_definition( 1, name, dtype, arity );
    EXPECT_EQ( _T("Color"), name );
    EXPECT_EQ( frantic::channels::data_type_float16, dtype );
    EXPECT_EQ( 3u, arity );
    cm.get_channel_definition( 2, name, dtype, arity );
    EXPECT_EQ( _T("Density"), name );
    EXPECT_EQ( frantic::channels::data_type_float32, dtype );
    EXPECT_EQ( 1u, arity );

    // The global metadata
    const property_map& gm = prt2.get_general_metadata();
    EXPECT_TRUE( gm.has_property( _T("LengthUnitInMicrometers") ) );
    EXPECT_EQ( 25400.0, gm.get<double>( _T("LengthUnitInMicrometers") ) );
    EXPECT_TRUE( gm.has_property( _T("CoordSys") ) );
    EXPECT_EQ( frantic::graphics::coordinate_system::right_handed_zup, gm.get<int32_t>( _T("CoordSys") ) );
    EXPECT_TRUE( gm.has_property( _T("FrameRate") ) );
    pair<uint32_t, uint32_t> frameRate = gm.get<pair<uint32_t, uint32_t>>( _T("FrameRate") );
    EXPECT_EQ( 30000U, frameRate.first );
    EXPECT_EQ( 1001U, frameRate.second );
    const channel_map& gmcm = gm.get_channel_map();
    EXPECT_EQ( 3, gmcm.channel_count() );
    gmcm.get_channel_definition( _T("LengthUnitInMicrometers"), dtype, arity );
    EXPECT_EQ( frantic::channels::data_type_float64, dtype );
    EXPECT_EQ( 1u, arity );
    gmcm.get_channel_definition( _T("CoordSys"), dtype, arity );
    EXPECT_EQ( frantic::channels::data_type_int32, dtype );
    EXPECT_EQ( 1u, arity );
    gmcm.get_channel_definition( _T("FrameRate"), dtype, arity );
    EXPECT_EQ( frantic::channels::data_type_uint32, dtype );
    EXPECT_EQ( 2u, arity );

    // The Position channel metadata
    const property_map& pcm = prt2.get_channel_metadata( _T("Position") );
    EXPECT_TRUE( pcm.has_property( _T("Extents") ) );
    EXPECT_EQ( boundbox3f( -1, 1, -1, 1, -1, 1 ), pcm.get<boundbox3f>( _T("Extents") ) );
    EXPECT_TRUE( pcm.has_property( _T("Interpretation") ) );
    EXPECT_EQ( _T("Point"), pcm.get_cvt<frantic::tstring>( _T("Interpretation") ) );
    const channel_map& pcmcm = pcm.get_channel_map();
    EXPECT_EQ( 2u, pcmcm.channel_count() );
    pcmcm.get_channel_definition( _T("Extents"), dtype, arity );
    EXPECT_EQ( frantic::channels::data_type_float32, dtype );
    EXPECT_EQ( 6u, arity );
    pcmcm.get_channel_definition( _T("Interpretation"), dtype, arity );
    EXPECT_EQ( frantic::channels::data_type_string, dtype );
    EXPECT_EQ( 1u, arity );

    // The combined metadata
    // The metadata pointers aren't required to be equal, but they are
    // in the current implementation.  I'm taking advantage of this so
    // I can compare the pointers instead of the individual data items.
    const frantic::particles::particle_file_metadata& metadata = prt2.get_metadata();
    EXPECT_EQ( &prt2.get_general_metadata(), &metadata.get_general_metadata() );
    EXPECT_EQ( &prt2.get_channel_metadata( _T("Position") ), metadata.get_channel_metadata( _T("Position") ) );

    // No other channel metadata
    EXPECT_EQ( 0u, prt2.get_channel_metadata( _T("Color") ).get_channel_map().channel_count() );
    EXPECT_EQ( 0u, prt2.get_channel_metadata( _T("Density") ).get_channel_map().channel_count() );

    // Verify the particle chunk index
    pair<int64_t, int64_t> ex;
    ex = prt2.get_particle_chunk_extents( 0 );
    EXPECT_EQ( 0, ex.first );
    EXPECT_EQ( 2, ex.second );
    ex = prt2.get_particle_chunk_extents( 1 );
    EXPECT_EQ( 2, ex.first );
    EXPECT_EQ( 8, ex.second );

    // Check the contents of the particles themselves
    particle_array parr( prt2.get_channel_map() );
    channel_cvt_accessor<vector3f> posAccessor = prt2.get_channel_map().get_cvt_accessor<vector3f>( _T("Position") );
    prt2.read_particle_chunk( 0, parr );
    EXPECT_EQ( vector3f( -1, -1, -1 ), posAccessor( parr[0] ) );
    EXPECT_EQ( vector3f( -1, -1, 1 ), posAccessor( parr[1] ) );
    EXPECT_EQ( 2u, parr.particle_count() );
    prt2.read_particle_chunk( 1, parr );
    EXPECT_EQ( 6u, parr.particle_count() );
    EXPECT_EQ( vector3f( -1, 1, -1 ), posAccessor( parr[0] ) );
    EXPECT_EQ( vector3f( -1, 1, 1 ), posAccessor( parr[1] ) );
    EXPECT_EQ( vector3f( 1, -1, -1 ), posAccessor( parr[2] ) );
    EXPECT_EQ( vector3f( 1, -1, 1 ), posAccessor( parr[3] ) );
    EXPECT_EQ( vector3f( 1, 1, -1 ), posAccessor( parr[4] ) );
    EXPECT_EQ( vector3f( 1, 1, 1 ), posAccessor( parr[5] ) );
}

TEST( PRT2, WriteHandcraftedFile ) {
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
    prt2.write_general_metadata_filechunk<double>( _T("LengthUnitInMicrometers"), 25400.0 );
    prt2.write_general_metadata_filechunk<int32_t>( _T("CoordSys"), 2 );
    prt2.write_general_metadata_filechunk<pair<uint32_t, uint32_t>>( _T("FrameRate"), make_pair( 30000, 1001 ) );
    prt2.write_channel_metadata_filechunk<frantic::tstring>( _T("Position"), _T("Interpretation"), _T("Point") );

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
                chunk->particleCount = 2;
                chunk->uncompressed.resize( m_particles.get_channel_map().structure_size() * 2 );
                memcpy( &chunk->uncompressed[0], m_particles.at( 0 ), chunk->uncompressed.size() );
                return chunk;
            }
            case 1: {
                prt2_writer::particle_chunk* chunk = new prt2_writer::particle_chunk();
                chunk->particleCount = 6;
                chunk->uncompressed.resize( m_particles.get_channel_map().structure_size() * 6 );
                memcpy( &chunk->uncompressed[0], m_particles.at( 2 ), chunk->uncompressed.size() );
                return chunk;
            }
            case 2: {
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
    prt2.write_particle_chunks( chunkGen, 8, nullProgress, false, frantic::prtfile::prt2_compression_uncompressed );

    // TODO: Add ability of the prt2_writer to automatically produce these extents
    prt2.write_channel_metadata_filechunk<boundbox3f>( _T("Position"), _T("Extents"),
                                                       boundbox3f( -1, 1, -1, 1, -1, 1 ) );

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

TEST( PRT2, CleanupTempFileOnError ) {
    using namespace frantic::files;

    channel_map pcm;
    pcm.end_channel_definition();
    pcm.append_channel<vector3f>( _T("Position") );

    boost::filesystem::path tempPath = boost::filesystem::unique_path(
        boost::filesystem::temp_directory_path() / boost::filesystem::path( "%%%%-%%%%-%%%%-%%%%" ) );
    boost::filesystem::create_directory( tempPath );

    scoped_file_cleanup cleanup;
    cleanup.add( tempPath );

    boost::filesystem::path targetFile =
        boost::filesystem::unique_path( tempPath / boost::filesystem::path( "%%%%-%%%%-%%%%-%%%%.prt" ) );

    try {
        prt2_writer prt2;
        prt2.open( frantic::files::to_tstring( targetFile ), pcm, true, tempPath );
        throw std::runtime_error( "drop due to exception." );
    } catch( std::exception& /*e*/ ) {
        // pass
    }

    std::vector<boost::filesystem::path> paths;
    paths.assign( boost::filesystem::directory_iterator( tempPath ), boost::filesystem::directory_iterator() );

    EXPECT_EQ( 0, paths.size() );
}
