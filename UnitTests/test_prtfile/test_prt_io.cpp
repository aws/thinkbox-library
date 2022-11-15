// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

// Problem with C++11 and linking Boost Filesystem -> undefined reference to copy_file(). This function is using
// BOOST_SCOPED_ENUMS and the fix is to add the following line.
#ifdef __linux__
#define BOOST_NO_CXX11_SCOPED_ENUMS
#endif

#include <frantic/particles/streams/prt_particle_istream.hpp>
#include <frantic/particles/streams/prt_particle_ostream.hpp>

#include <frantic/files/paths.hpp>
#include <frantic/files/scoped_file_cleanup.hpp>

#include <frantic/channels/channel_map.hpp>

#include <frantic/graphics/vector3f.hpp>

#include <frantic/particles/particle_file_stream_factory.hpp>

using frantic::files::scoped_file_cleanup;
namespace fs = boost::filesystem;

TEST( PRT, CleanupTempFileOnError ) {
    using namespace frantic::files;
    using namespace frantic::channels;
    using namespace frantic::graphics;
    using namespace frantic::particles;

    channel_map pcm;
    pcm.end_channel_definition();
    pcm.append_channel<vector3f>( _T("Position") );

    fs::path tempPath = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%" );
    fs::create_directory( tempPath );

    scoped_file_cleanup cleanup;
    cleanup.add( tempPath );

    fs::path targetFile = fs::unique_path( tempPath / fs::path( "%%%%-%%%%-%%%%-%%%%.prt" ) );

    try {
        streams::prt_particle_ostream outputStream( frantic::files::to_tstring( targetFile ), pcm, pcm, -1, -1,
                                                    tempPath );
        throw std::runtime_error( "drop due to exception." );
    } catch( std::exception& /*e*/ ) {
        // pass
    }

    std::vector<fs::path> paths;
    paths.assign( fs::directory_iterator( tempPath ), fs::directory_iterator() );

    EXPECT_EQ( 0, paths.size() );
}

TEST( PRT, ParticleStreamHeader ) {
    fs::path path = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.prt" );
    scoped_file_cleanup fileCleanup;
    fileCleanup.add( path );

    frantic::channels::channel_map channelMap;
    channelMap.define_channel<frantic::graphics::vector3f>( _T("Position") );
    channelMap.end_channel_definition();

    {
        frantic::particles::particle_ostream_ptr pout( new frantic::particles::streams::prt_particle_ostream(
            frantic::files::to_tstring( path ), channelMap, channelMap ) );

        std::vector<char> buffer( channelMap.structure_size() );

        frantic::channels::channel_accessor<frantic::graphics::vector3f> acc(
            channelMap.get_accessor<frantic::graphics::vector3f>( _T("Position") ) );

        acc( buffer ).set( 0, 0, 0 );
        pout->put_particle( buffer );

        acc( buffer ).set( 1, 1, 1 );
        pout->put_particle( buffer );
    }

    {
        boost::shared_ptr<frantic::particles::streams::prt_particle_istream> pin(
            new frantic::particles::streams::prt_particle_istream( frantic::files::to_tstring( path ) ) );

        EXPECT_EQ( pin->particle_count(), 2 );

        const frantic::channels::property_map* metadataPtr = pin->get_channel_metadata( _T("Position") );
        ASSERT_TRUE( metadataPtr != 0 );

        const frantic::channels::property_map& metadata = *metadataPtr;
        ASSERT_TRUE( metadata.has_property( _T("Extents") ) );

        const frantic::graphics::boundbox3f actualBox = metadata.get<frantic::graphics::boundbox3f>( _T("Extents") );
        const frantic::graphics::boundbox3f expectedBox( frantic::graphics::vector3f( 0 ),
                                                         frantic::graphics::vector3f( 1 ) );

        EXPECT_EQ( expectedBox.minimum(), actualBox.minimum() );
        EXPECT_EQ( expectedBox.maximum(), actualBox.maximum() );
    }
}

// Non-standard metadata with property name "Position.Extents" should be
// interpreted as a per-channel property for the "Position" channel named
// "Extents".
TEST( PRT, ChannelDotPropertyMetadataWorkaround ) {
    const fs::path filename = fs::unique_path( fs::temp_directory_path() / "%%%%-%%%%-%%%%-%%%%.prt" );

    frantic::files::scoped_file_cleanup cleanup;
    cleanup.add( filename );

    {
        // clang-format off
		static unsigned const char chunkSection[] = {
			// chunk type
			'M', 'e', 't', 'a',
			// data length
			46, 0, 0, 0,
			// channel name
			0,
			// name
			'P', 'o', 's', 'i', 't', 'i', 'o', 'n', '.', 'E', 'x', 't', 'e', 'n', 't', 's', 0,
			// type (float32)
			4, 0, 0, 0,
			// data
			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
			'S', 't', 'o', 'p',
			0, 0, 0, 0 };

		static unsigned const char header[] = {
			192, 'P', 'R', 'T', '\r', '\n', 26, '\n',
			56 + sizeof( chunkSection ), 0, 0, 0,
			'E', 'x', 't', 'e', 'n', 's', 'i', 'b', 'l', 'e', ' ', 'P', 'a', 'r', 't', 'i', 'c', 'l', 'e', ' ', 'F', 'o', 'r', 'm', 'a', 't', 0, 0, 0, 0, 0, 0,
			2, 0, 0, 0,
			1, 0, 0, 0, 0, 0, 0, 0 };

		static unsigned const char reservedSection[] = {
			4, 0, 0, 0 };

		static unsigned const char channelSection[] = {
			// number of channels
			1, 0, 0, 0,
			44, 0, 0, 0,
			'P', 'o', 's', 'i', 't', 'i', 'o', 'n', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			// channel type (float32)
			4, 0, 0, 0,
			// arity
			3, 0, 0, 0,
			// offset
			0, 0, 0, 0 };
        // clang-format on

        std::ofstream out( filename.c_str(), std::ios::out | std::ios::binary );

        out.write( reinterpret_cast<const char*>( header ), sizeof( header ) );
        out.write( reinterpret_cast<const char*>( chunkSection ), sizeof( chunkSection ) );
        out.write( reinterpret_cast<const char*>( reservedSection ), sizeof( reservedSection ) );
        out.write( reinterpret_cast<const char*>( channelSection ), sizeof( channelSection ) );

        {
            frantic::files::zlib_deflate_ostream outDeflate( out, _T("memory") );

            const char particle[12] = { 0 };

            outDeflate.write( particle, sizeof( particle ) );
        }
    }

    frantic::particles::particle_file_stream_factory_object factory;
    frantic::particles::particle_file_metadata metadata;
    frantic::particles::particle_istream_ptr pin =
        factory.create_istream( frantic::files::to_tstring( filename ), &metadata );

    frantic::channels::property_map* positionMetadata = metadata.get_channel_metadata( _T("Position") );
    ASSERT_TRUE( positionMetadata != NULL );
    EXPECT_TRUE( positionMetadata->has_property( _T("Extents") ) );
}
