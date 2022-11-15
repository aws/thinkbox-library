// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#ifdef _WIN32
#include <frantic/win32/utility.hpp>
#include <windows.h>
#endif

#include <frantic/logging/logging_level.hpp>

#include <frantic/files/compression_stream.hpp>
#include <frantic/files/paths.hpp>

#include <frantic/volumetrics/levelset/rle_level_set_file_io.hpp>

#include <boost/numeric/conversion/bounds.hpp>
#include <frantic/files/files.hpp>

#include <frantic/graphics/vector3f.hpp>

using frantic::graphics::vector3f;

using namespace std;
using namespace boost;
using namespace frantic;
using namespace frantic::files;

namespace frantic {
namespace volumetrics {
namespace levelset {

// This general function loads any type of level set input file.  It should be extended when new
// file format support is added.

void read_rle_level_set_file( const frantic::tstring& file, rle_level_set& outLevelSet ) {
    frantic::tstring extension = strings::to_lower( files::extension_from_path( file ) );
    if( extension == _T(".rls") ) {
        read_rls_rle_level_set_file( file, outLevelSet );
    } else {
        throw runtime_error( "read_rle_level_set_file() - The file extension of the filename given, \"" +
                             frantic::strings::to_string( file ) + "\", is not a recognized level set file type." );
    }
}

void read_rle_level_set_file_header( const frantic::tstring& file, rle_level_set& outLevelSet ) {
    frantic::tstring extension = strings::to_lower( files::extension_from_path( file ) );
    if( extension == _T(".rls") ) {
        read_rls_rle_level_set_file_header( file, outLevelSet );
    } else {
        throw runtime_error( "read_rle_level_set_file_header() - The file extension of the filename given, \"" +
                             frantic::strings::to_string( file ) + "\", is not a recognized level set file type." );
    }
}

namespace detail {
#pragma pack( 4 )
struct rls_header {
    char magicNumber[8];
    char fileFormatID[32];
    boost::int32_t versionNumber;
    float worldOriginX, worldOriginY, worldOriginZ;
    float voxelLength;
    float interfaceVoxelWidthInside, interfaceVoxelWidthOutside;
};

struct rls_channel_header {
    char channelName[32];
    char channelDataType[16];
    boost::int32_t arity;
    boost::int64_t channelDataSeekLocation;
};

struct ris_header {
    boost::int32_t xmin, xmax, ymin, ymax, zmin, zmax;
    boost::int32_t exteriorRegionCode;
    boost::uint32_t bcToRunIndexSize, runIndexDataSize, dataSize;
};

#pragma pack()

void read_header( const frantic::tstring& file, ifstream& fin, rls_header& fileHeader,
                  vector<rls_channel_header>& channelHeaders ) {

    ////////
    // Read and validate the header
    ////////

    fin.read( reinterpret_cast<char*>( &fileHeader ), sizeof( detail::rls_header ) );
    if( !fin )
        throw runtime_error( "read_rls_rle_level_set_file() - Error reading the header of file \"" +
                             frantic::strings::to_string( file ) + "\"." );

    static unsigned char rlsMagicNumber[8] = { 192, 'R', 'L', 'S', '\r', '\n', 26, '\n' };
    if( memcmp( rlsMagicNumber, fileHeader.magicNumber, 8 ) != 0 )
        throw runtime_error( "read_rls_rle_level_set_file() - The magic number of file \"" +
                             frantic::strings::to_string( file ) + "\" was not correct." );

    if( strncmp( fileHeader.fileFormatID, "RLE Level Set File", 31 ) != 0 )
        throw runtime_error( "read_rls_rle_level_set_file() - The file format ID string of file \"" +
                             frantic::strings::to_string( file ) + "\" was not correct." );

    if( fileHeader.versionNumber != 1 )
        throw runtime_error(
            "read_rls_rle_level_set_file() - The version number \"" + frantic::strings::to_string( file ) + "\" was " +
            lexical_cast<string>( fileHeader.versionNumber ) + ", but only version 1 is supported in this program." );

    ////////
    // Read the channel header
    ////////

    boost::int32_t channelCount;
    fin.read( reinterpret_cast<char*>( &channelCount ), 4 );
    if( !fin )
        throw runtime_error( "read_rls_rle_level_set_file() - Error reading the channel header of file \"" +
                             frantic::strings::to_string( file ) + "\"." );

    if( channelCount < 1 || channelCount > 10000 )
        throw runtime_error( "read_rls_rle_level_set_file() - The file \"" + frantic::strings::to_string( file ) +
                             "\" specified that there are " + lexical_cast<string>( channelCount ) +
                             " channels, which seems like an unreasonable value to this file reading code." );

    channelHeaders.resize( channelCount );
    fin.read( reinterpret_cast<char*>( &channelHeaders[0] ), sizeof( detail::rls_channel_header ) * channelCount );
    if( !fin )
        throw runtime_error( "read_rls_rle_level_set_file() - Error reading the channel header of file \"" +
                             frantic::strings::to_string( file ) + "\"." );
}

void read_rle_index_spec( const frantic::tstring& file, ifstream& fin, rle_index_spec& outRis ) {
    ////////
    // Read the rle index spec
    ////////

    zlib_inflate_istream zlibIn;
    zlibIn.open( fin, file );

    detail::ris_header risHeader;
    zlibIn.read( reinterpret_cast<char*>( &risHeader ), sizeof( detail::ris_header ) );

    vector<boost::int32_t> bcToRunIndex( risHeader.bcToRunIndexSize );
    vector<run_data> runIndexData( risHeader.runIndexDataSize );
    if( risHeader.bcToRunIndexSize != 0 )
        zlibIn.read( reinterpret_cast<char*>( &bcToRunIndex[0] ),
                     sizeof( boost::int32_t ) * risHeader.bcToRunIndexSize );
    if( risHeader.runIndexDataSize != 0 )
        zlibIn.read( reinterpret_cast<char*>( &runIndexData[0] ), sizeof( run_data ) * risHeader.runIndexDataSize );

    zlibIn.close();

    frantic::graphics::boundbox3 outerBounds( risHeader.xmin, risHeader.xmax, risHeader.ymin, risHeader.ymax,
                                              risHeader.zmin, risHeader.zmax );
    // We always want to validate consistency when loading from a file.
    outRis.set_with_swap( /*checkConsistency=*/true, outerBounds, risHeader.dataSize, risHeader.exteriorRegionCode,
                          bcToRunIndex, runIndexData, file );
}

template <class RLEObject>
void read_rle_named_channels( const frantic::tstring& file, ifstream& fin,
                              const channels::channel_propagation_policy& policy,
                              const vector<detail::rls_channel_header>& channelHeaders, RLEObject& result ) {
    zlib_inflate_istream zlibIn;

    for( unsigned i = 0; i < channelHeaders.size(); ++i ) {

        // see if this is in the include
        //			bool includeThisChannel = false;

        // check if this channel shound be included
        if( policy.is_channel_included( frantic::strings::to_tstring( channelHeaders[i].channelName ) ) ) {
            // Get the data type
            channels::data_type_t dataType =
                channels::channel_data_type_from_string( channelHeaders[i].channelDataType );

            // Validate the channel information
            if( channelHeaders[i].channelDataSeekLocation == 0 )
                throw runtime_error( string() + "read_rle_named_channels() - The seek location for channel \"" +
                                     channelHeaders[i].channelName + "\" in the file \"" +
                                     frantic::strings::to_string( file ) +
                                     "\" was 0, indicating in incompletely written file." );
            if( !channels::is_valid_channel_name( channelHeaders[i].channelName ) )
                throw runtime_error(
                    string() + "read_rle_named_channels() - The channel name \"" + channelHeaders[i].channelName +
                    "\" in the file \"" + frantic::strings::to_string( file ) +
                    "\" is not valid.  Channel names must start with a letter or underscore, and consist "
                    "entirey of letters, numbers, and underscores." );
            if( dataType == channels::data_type_invalid )
                throw runtime_error( string() + "read_rle_named_channels() - The channel with the name \"" +
                                     channelHeaders[i].channelName + "\" in the file \"" +
                                     frantic::strings::to_string( file ) + "\" has an invalid data type, \"" +
                                     channelHeaders[i].channelDataType + "\"." );
            if( channelHeaders[i].arity < 1 )
                throw runtime_error( string() + "read_rle_named_channels() - The channel with the name \"" +
                                     channelHeaders[i].channelName + "\" in the file \"" +
                                     frantic::strings::to_string( file ) + "\" has an invalid arity " +
                                     lexical_cast<string>( channelHeaders[i].arity ) +
                                     ".  Its arity must be positive." );

            // Check whether the channel already exists in the file
            if( result.has_channel( frantic::strings::to_tstring( channelHeaders[i].channelName ) ) )
                throw runtime_error(
                    string() + "read_rle_named_channels() - The channel \"" + channelHeaders[i].channelName +
                    "\" appears more than once in the channel header for file \"" +
                    frantic::strings::to_string( file ) + "\", all channels should have unique names." );

            result.add_channel( frantic::strings::to_tstring( channelHeaders[i].channelName ), channelHeaders[i].arity,
                                dataType );
            rle_channel_general_accessor channelAccessor =
                result.get_channel_general_accessor( frantic::strings::to_tstring( channelHeaders[i].channelName ) );

            if( result.get_rle_index_spec().data_size() > 0 ) {
                // Seek to the location within the file
                fin.clear();
                fin.seekg( (std::streamoff)channelHeaders[i].channelDataSeekLocation );
                if( !fin )
                    throw runtime_error( "read_rle_named_channels() - Failed to seek in file \"" +
                                         frantic::strings::to_string( file ) + "\" to location " +
                                         lexical_cast<string>( channelHeaders[i].channelDataSeekLocation ) +
                                         " in order to read the \"" + channelHeaders[i].channelName + "\" channel." );

                // Read in the data from the compressed stream
                zlibIn.open( fin, file );
                FF_LOG( debug ) << "Opened Zlib Stream for Channels" << endl;
                zlibIn.read( channelAccessor.data( 0 ),
                             int( channelAccessor.primitive_size() * result.get_rle_index_spec().data_size() ) );
                FF_LOG( debug ) << "read Channel Data" << endl;
                zlibIn.close();
            }
        }
    }
}
} // namespace detail

void read_rls_rle_level_set_file_header( const frantic::tstring& file, rle_level_set& outLevelSet ) {

    ifstream fin( file.c_str(), ios::in | ios::binary );

    if( !fin.is_open() )
        throw runtime_error( "read_rls_rle_level_set_file() - Failed to open the file \"" +
                             frantic::strings::to_string( file ) + "\" for input." );

    detail::rls_header fileHeader;
    vector<detail::rls_channel_header> channelHeaders;

    detail::read_header( file, fin, fileHeader, channelHeaders );

    // Create an empty level set with the read in header info
    vector3f origin( fileHeader.worldOriginX, fileHeader.worldOriginY, fileHeader.worldOriginZ );
    rle_index_spec ris;
    vector<float> distanceData;
    voxel_coord_system vcs( origin, fileHeader.voxelLength );
    outLevelSet.set_with_swap( vcs, ris, distanceData, fileHeader.interfaceVoxelWidthInside,
                               fileHeader.interfaceVoxelWidthOutside );

    // Add the channels
    for( size_t i = 0; i < channelHeaders.size(); i++ )
        outLevelSet.add_channel( frantic::strings::to_tstring( channelHeaders[i].channelName ), channelHeaders[i].arity,
                                 channels::channel_data_type_from_string( channelHeaders[i].channelDataType ) );
}

void read_rls_rle_level_set_file( const frantic::tstring& file, rle_level_set& outLevelSet ) {
    ifstream fin( file.c_str(), ios::in | ios::binary );
    zlib_inflate_istream zlibIn;

    if( !fin.is_open() )
        throw runtime_error( "read_rls_rle_level_set_file() - Failed to open the file \"" +
                             frantic::strings::to_string( file ) + "\" for input." );

    detail::rls_header fileHeader;
    vector<detail::rls_channel_header> channelHeaders;

    detail::read_header( file, fin, fileHeader, channelHeaders );

    ////////
    // Read the rle index spec
    ////////

    rle_index_spec ris;
    detail::read_rle_index_spec( file, fin, ris );

    ////////
    // Read the "SignedDistance" channel
    ////////

    detail::rls_channel_header* signedDistanceChannelHeader = 0;
    // Find the "SignedDistance" channel in the headers
    for( unsigned i = 0; i < channelHeaders.size(); ++i ) {
        if( strncmp( channelHeaders[i].channelName, "SignedDistance", 31 ) == 0 ) {
            if( signedDistanceChannelHeader == 0 )
                signedDistanceChannelHeader = &channelHeaders[i];
            else
                throw runtime_error(
                    "read_rls_rle_level_set_file() - The channel \"SignedDistance\" appears more than once in "
                    "the channel header for file \"" +
                    frantic::strings::to_string( file ) + "\", all channels should have unique names." );
        }
    }

    if( signedDistanceChannelHeader == 0 )
        throw runtime_error(
            "read_rls_rle_level_set_file() - There is no channel with the name \"SignedDistance\" in the file \"" +
            frantic::strings::to_string( file ) + "\"." );

    if( signedDistanceChannelHeader->channelDataSeekLocation == 0 )
        throw runtime_error(
            "read_rls_rle_level_set_file() - The seek location for channel \"SignedDistance\" in the file \"" +
            frantic::strings::to_string( file ) + "\" was 0, indicating in incompletely written file." );

    channels::data_type_t signedDistanceDataType =
        channels::channel_data_type_from_string( signedDistanceChannelHeader->channelDataType );

    if( signedDistanceDataType != channels::data_type_float32 )
        throw runtime_error(
            "read_rls_rle_level_set_file() - The channel with the name \"SignedDistance\" in the file \"" +
            frantic::strings::to_string( file ) +
            "\" has the incorrect data type.  Its data type should be \"float32\"." );
    if( signedDistanceChannelHeader->arity != 1 )
        throw runtime_error(
            "read_rls_rle_level_set_file() - The channel with the name \"SignedDistance\" in the file \"" +
            frantic::strings::to_string( file ) + "\" has the incorrect arity.  Its arity should be 1." );

    // Now we've got the channel, and confirmed that it is a float32 channel with arity 1, so we can read it directly
    // into a distanceData array.

    // TODO: In the future we will want to add support for converting from float16 or float64.

    vector<float> distanceData( ris.data_size() );

    if( ris.data_size() > 0 ) {
        // Seek to the location within the file
        fin.clear();
        fin.seekg( (std::streamoff)signedDistanceChannelHeader->channelDataSeekLocation );
        if( !fin )
            throw runtime_error( "read_rls_rle_level_set_file() - Failed to seek in file \"" +
                                 frantic::strings::to_string( file ) + "\" to location " +
                                 lexical_cast<string>( signedDistanceChannelHeader->channelDataSeekLocation ) +
                                 " in order to read the \"SignedDistance\" channel." );

        // Read in the data from the compressed stream
        zlibIn.open( fin, file );
        zlibIn.read( reinterpret_cast<char*>( &distanceData[0] ), int( sizeof( float ) * ris.data_size() ) );
        zlibIn.close();
    }

    ////////
    // Build the basic rle level set
    ////////

    voxel_coord_system vcs( vector3f( fileHeader.worldOriginX, fileHeader.worldOriginY, fileHeader.worldOriginZ ),
                            fileHeader.voxelLength );

    // The level set is read into this temporary variable first, then swapped into outLevelSet.  This is
    // so that if an exception is thrown during reading, the outLevelSet is left untouched.
    rle_level_set result;
    result.set_with_swap( vcs, ris, distanceData, fileHeader.interfaceVoxelWidthInside,
                          fileHeader.interfaceVoxelWidthOutside );

    ////////
    // Fill in the rest of the channels
    ////////

    // exclude the signed distance channel from propagating
    // since this object is an rle_level_set and has already loaded it
    channels::channel_propagation_policy policy;
    policy.set_to_exclude_policy();
    policy.add_channel( _T("SignedDistance") );

    detail::read_rle_named_channels( file, fin, policy, channelHeaders, result );

    ////////
    // Swap the result into outLevelSet
    ////////
    result.swap( outLevelSet );
}

void read_rle_voxel_field_file( const frantic::tstring& file, fluids::rle_voxel_field& outVoxelField ) {
    // call the main function with a empty, fully inclusive channel policy
    read_rle_voxel_field_file( file, channels::channel_propagation_policy(), outVoxelField );
}

void read_rle_voxel_field_file( const frantic::tstring& file, const channels::channel_propagation_policy& policy,
                                fluids::rle_voxel_field& outVoxelField ) {
    ifstream fin( file.c_str(), ios::in | ios::binary );
    zlib_inflate_istream zlibIn;

    if( !fin.is_open() )
        throw runtime_error( "read_rls_rle_level_set_file() - Failed to open the file \"" +
                             frantic::strings::to_string( file ) + "\" for input." );

    detail::rls_header fileHeader;
    vector<detail::rls_channel_header> channelHeaders;

    detail::read_header( file, fin, fileHeader, channelHeaders );

    ////////
    // Read the rle index spec
    ////////

    rle_index_spec ris;
    detail::read_rle_index_spec( file, fin, ris );

    ////////
    // Build the basics for the rle voxel field
    ////////

    voxel_coord_system vcs( vector3f( fileHeader.worldOriginX, fileHeader.worldOriginY, fileHeader.worldOriginZ ),
                            fileHeader.voxelLength );

    // The level ste is read into this temporary variable first, then swapped into outLevelSet.  This is
    // so that if an exception is thrown during reading, the outLevelSet is left untouched.
    fluids::rle_voxel_field result;
    result.set_with_swap( vcs, ris );

    ////////
    // Fill in the rest of the channels
    ////////
    detail::read_rle_named_channels( file, fin, policy, channelHeaders, result );

    ////////
    // Swap the result into outLevelSet
    ////////
    result.swap( outVoxelField );
}

namespace detail {
class rls_writer_tempfile_deleter {
    frantic::tstring filename;

  public:
    rls_writer_tempfile_deleter( const frantic::tstring& filenameToDelete )
        : filename( filenameToDelete ) {}

    ~rls_writer_tempfile_deleter() {
        // DeleteFile( filename.c_str() );
        frantic::files::delete_file( filename.c_str() );
    }
};

/**
 * Writes an rle header for rle containers to an output stream. The writeSignedDistanceFirst flag specifies whether
 * to write the 'SignedDistance' channel as the first channel. A requirement for an rle_level_set
 *
 *@param fout the output file stream
 *@param file the output filename
 *@param object the rle object (rle_voxel_field, rle_level_set)
 *@param writeSignedDistanceFirst whether or not to write the "SignedDistance" channel first. Required for rle_level_set
 */
template <class RLEObject>
void write_rle_header( ostream& fout, const frantic::tstring& file, const RLEObject& object,
                       const channels::channel_propagation_policy& policy, bool writeSignedDistanceFirst,
                       float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside ) {
    ////////
    // Write the main file header
    ////////

    detail::rls_header rlsHeader;
    memset( &rlsHeader, 0, sizeof( detail::rls_header ) );

    // Fill in the rls_header structure
    static unsigned char rlsMagicNumber[8] = { 192, 'R', 'L', 'S', '\r', '\n', 26, '\n' };
    memcpy( rlsHeader.magicNumber, rlsMagicNumber, 8 );
    strncpy( rlsHeader.fileFormatID, "RLE Level Set File", 31 );
    rlsHeader.versionNumber = 1;
    rlsHeader.worldOriginX = object.get_voxel_coord_system().world_origin().x;
    rlsHeader.worldOriginY = object.get_voxel_coord_system().world_origin().y;
    rlsHeader.worldOriginZ = object.get_voxel_coord_system().world_origin().z;
    rlsHeader.voxelLength = object.get_voxel_coord_system().voxel_length();
    rlsHeader.interfaceVoxelWidthInside = interfaceVoxelWidthInside;
    rlsHeader.interfaceVoxelWidthOutside = interfaceVoxelWidthOutside;

    ////////////
    //// Write the header
    //////////
    fout.write( reinterpret_cast<const char*>( &rlsHeader ), sizeof( detail::rls_header ) );
    if( !fout )
        throw std::runtime_error( "write_rle_header() - Failed to write the RLS file header to file \"" +
                                  frantic::strings::to_string( file ) + "\"." );

    //////////
    //// Write the Channels Header
    //////////

    //// Get all the channel names
    vector<frantic::tstring> channelNames;
    object.get_channel_names( channelNames );
    policy.filter_channel_vector( channelNames ); // apply channel propagation policy

    boost::int32_t channelCount = boost::int32_t( channelNames.size() );

    // the rls level set that uses "writeSignedDistanceFirst" will always count have a signed distance channel regarless
    // of propagation policy
    if( writeSignedDistanceFirst ) {
        ++channelCount;
    }

    fout.write( reinterpret_cast<const char*>( &channelCount ), 4 );

    if( !fout )
        throw runtime_error( "write_rle_header() - Error writing the channel header to file \"" +
                             frantic::strings::to_string( file ) + "\"." );

    detail::rls_channel_header rlsChannelHeader;

    // the rls level set that uses "writeSignedDistanceFirst" needs to always write a signed distance regardless of
    // propagation policy
    if( writeSignedDistanceFirst ) {
        // Write the "SignedDistance" channel as channel 0
        memset( &rlsChannelHeader, 0, sizeof( detail::rls_channel_header ) );
        strncpy( rlsChannelHeader.channelName, "SignedDistance", 31 );
        strncpy( rlsChannelHeader.channelDataType,
                 frantic::strings::to_string( channels::channel_data_type_str( channels::data_type_float32 ) ).c_str(),
                 15 );
        rlsChannelHeader.arity = 1;
        fout.write( reinterpret_cast<const char*>( &rlsChannelHeader ), sizeof( detail::rls_channel_header ) );
        if( !fout )
            throw runtime_error( "write_rle_header() - Error writing the \"SignedDistance\" channel header to file \"" +
                                 frantic::strings::to_string( file ) + "\"." );
    }

    // Loop through all the named channels and write out header entries for them
    for( unsigned i = 0; i < channelNames.size(); ++i ) {
        const_rle_channel_general_accessor channelAccessor = object.get_channel_general_accessor( channelNames[i] );

        memset( &rlsChannelHeader, 0, sizeof( detail::rls_channel_header ) );
        strncpy( rlsChannelHeader.channelName, frantic::strings::to_string( channelNames[i] ).c_str(), 31 );
        strncpy( rlsChannelHeader.channelDataType,
                 frantic::strings::to_string( channels::channel_data_type_str( channelAccessor.data_type() ) ).c_str(),
                 15 );
        rlsChannelHeader.arity = (boost::int32_t)channelAccessor.arity();
        fout.write( reinterpret_cast<const char*>( &rlsChannelHeader ), sizeof( detail::rls_channel_header ) );
        if( !fout )
            throw runtime_error( "write_rle_header() - Error writing the \"" +
                                 frantic::strings::to_string( channelNames[i] ) + "\" channel header to file \"" +
                                 frantic::strings::to_string( file ) + "\"." );
    }
}

/**
 * Writes an rle index spec to an output stream.
 *
 *@param fout the output file stream
 *@param file the output filename
 *@param ris the rle index spec
 */
void write_rle_index_spec( ofstream& fout, const frantic::tstring& file, const rle_index_spec& ris ) {
    zlib_deflate_ostream zlibOut;
    zlibOut.open( fout, file );

    // Get the rle index spec value reference, and const references to its big arrays
    const vector<boost::int32_t>& bcToRunIndex = ris.get_bc_to_run_index_vector();
    const vector<run_data>& runIndexData = ris.get_run_index_data_vector();

    // Fill in the header data
    detail::ris_header risHeader;
    memset( &risHeader, 0, sizeof( detail::ris_header ) );
    frantic::graphics::boundbox3 outerBounds = ris.outer_bounds();
    risHeader.xmin = outerBounds.xminimum();
    risHeader.xmax = outerBounds.xmaximum();
    risHeader.ymin = outerBounds.yminimum();
    risHeader.ymax = outerBounds.ymaximum();
    risHeader.zmin = outerBounds.zminimum();
    risHeader.zmax = outerBounds.zmaximum();
    risHeader.exteriorRegionCode = ris.get_exterior_region_code();
    risHeader.bcToRunIndexSize = (boost::uint32_t)bcToRunIndex.size();
    risHeader.runIndexDataSize = (boost::uint32_t)runIndexData.size();
    risHeader.dataSize = (boost::uint32_t)ris.data_size();

    // Write the header data
    zlibOut.write( reinterpret_cast<const char*>( &risHeader ), sizeof( detail::ris_header ) );

    // Write the BC to Run Index array
    if( !bcToRunIndex.empty() )
        zlibOut.write( reinterpret_cast<const char*>( &bcToRunIndex[0] ),
                       sizeof( boost::int32_t ) * risHeader.bcToRunIndexSize );

    // Write the Run Index Data array
    if( !runIndexData.empty() )
        zlibOut.write( reinterpret_cast<const char*>( &runIndexData[0] ),
                       sizeof( run_data ) * risHeader.runIndexDataSize );

    zlibOut.close();
}

/**
 * Writes the named channels of rle voxel field to an output file stream
 *
 *@param fout the output file stream
 *@param file the output filename
 *@param field the rle voxel field
 */
void write_rle_named_channels( ofstream& fout, const frantic::tstring& file, const fluids::rle_voxel_field& field,
                               const channels::channel_propagation_policy& policy ) {
    zlib_deflate_ostream zlibOut;

    vector<frantic::tstring> channelNames;
    field.get_channel_names( channelNames );
    policy.filter_channel_vector( channelNames ); // apply channel propagation policy

    boost::int32_t channelCount = boost::int32_t( channelNames.size() );

    const rle_index_spec& ris = field.get_rle_index_spec();

    // Save the channel data seek offsets in this array
    vector<boost::int64_t> channelDataSeekOffsets( channelCount );

    // Now go through the rest of the channels, and write their data
    for( unsigned i = 0; i < channelNames.size(); ++i ) {
        const_rle_channel_general_accessor channelAccessor = field.get_channel_general_accessor( channelNames[i] );

        channelDataSeekOffsets[i] = fout.tellp();
        if( ris.data_size() > 0 ) {
            zlibOut.open( fout, file );
            zlibOut.write( channelAccessor.data( 0 ), channelAccessor.primitive_size() * ris.data_size() );
            zlibOut.close();
        }
    }

    ////////
    // Patch up the seek offsets in the channel headers so the match up to the actual data we wrote
    ////////

    for( int i = 0; i < channelCount; ++i ) {
        fout.seekp( sizeof( detail::rls_header ) + 4 + i * sizeof( detail::rls_channel_header ) + 52 );
        fout.write( reinterpret_cast<const char*>( &channelDataSeekOffsets[i] ), sizeof( boost::int64_t ) );
        if( !fout )
            throw runtime_error( "write_rle_named_channels() - Error writing the channel data offsets to file \"" +
                                 frantic::strings::to_string( file ) + "\"." );
    }
}

/**
 * Writes the named channels of rle level set to an output file stream
 *
 *@param fout the output file stream
 *@param file the output filename
 *@param levelSet the rle level set
 */
void write_rle_named_channels( ofstream& fout, const frantic::tstring& file, const rle_level_set& levelSet ) {
    zlib_deflate_ostream zlibOut;

    vector<frantic::tstring> channelNames;
    levelSet.get_channel_names( channelNames );

    boost::int32_t channelCount = boost::int32_t( channelNames.size() + 1 );

    // Save the channel data seek offsets in this array
    vector<boost::int64_t> channelDataSeekOffsets( channelCount );

    const rle_index_spec& ris = levelSet.get_rle_index_spec();

    // First write the "SignedDistance" channel data
    channelDataSeekOffsets[0] = fout.tellp();
    if( ris.data_size() > 0 ) {
        zlibOut.open( fout, file );
        zlibOut.write( reinterpret_cast<const char*>( &( levelSet.get_distance_data() )[0] ),
                       sizeof( float ) * ris.data_size() );
        zlibOut.close();
    }

    // Now go through the rest of the channels, and write their data
    for( unsigned i = 0; i < channelNames.size(); ++i ) {
        const_rle_channel_general_accessor channelAccessor = levelSet.get_channel_general_accessor( channelNames[i] );

        channelDataSeekOffsets[i + 1] = fout.tellp();
        if( ris.data_size() > 0 ) {
            zlibOut.open( fout, file );
            zlibOut.write( channelAccessor.data( 0 ), channelAccessor.primitive_size() * ris.data_size() );
            zlibOut.close();
        }
    }

    ////////
    // Patch up the seek offsets in the channel headers so the match up to the actual data we wrote
    ////////

    for( int i = 0; i < channelCount; ++i ) {
        fout.seekp( sizeof( detail::rls_header ) + 4 + i * sizeof( detail::rls_channel_header ) + 52 );
        fout.write( reinterpret_cast<const char*>( &channelDataSeekOffsets[i] ), sizeof( boost::int64_t ) );
        if( !fout )
            throw runtime_error( "write_rle_named_channels() - Error writing the channel data offsets to file \"" +
                                 frantic::strings::to_string( file ) + "\"." );
    }
}

} // namespace detail

void write_rle_voxel_field_file( const frantic::tstring& file, const fluids::rle_voxel_field& voxelfield ) {
    // call main function with an empty, fully inclusive channel policy
    write_rle_voxel_field_file( file, channels::channel_propagation_policy(), voxelfield );
}

void write_rle_voxel_field_file( const frantic::tstring& file, const channels::channel_propagation_policy& policy,
                                 const fluids::rle_voxel_field& field ) {
    frantic::tstring tempFilename = file + _T(".partial");

    // This object will delete the temp file if an exception gets thrown.
    detail::rls_writer_tempfile_deleter tempFilenameDeleter( tempFilename );

    // Put the ofstream in a nested block
    {
        ofstream fout( tempFilename.c_str(), ios::out | ios::binary );

        if( !fout.is_open() )
            throw std::runtime_error( "write_rle_voxel_field_file() - Failed to open the temporary file \"" +
                                      frantic::strings::to_string( tempFilename ) + "\" for output." );

        ////////
        // Write the main file header
        ////////

        detail::write_rle_header( fout, file, field, policy, false, 5, 5 );

        ////////
        // Write the rle_index_spec
        ////////

        detail::write_rle_index_spec( fout, file, field.get_rle_index_spec() );

        ////////
        // Write all the channel data
        ////////

        detail::write_rle_named_channels( fout, file, field, policy );

        ////////
        // Close the file, and make sure that worked
        ////////

        fout.close();
        if( fout.bad() )
            throw std::runtime_error( "write_rle_voxel_field_file() - Failed to close the output stream for file \"" +
                                      frantic::strings::to_string( file ) + "\"." );
    }

    // Finally, Move the file we've been writing to its proper destination path
    // DeleteFile(file.c_str());
    // if( !MoveFile(tempFilename.c_str(), file.c_str()) )
    frantic::files::delete_file( file.c_str() );
    if( 0 != frantic::files::rename_file( tempFilename.c_str(), file.c_str() ) )
        throw std::runtime_error( "write_rle_voxel_field_file() - Failed to move the temporary file \"" +
                                  frantic::strings::to_string( tempFilename ) + "\" to destination \"" +
                                  frantic::strings::to_string( file ) + "\": "
#ifdef _WIN32
                                  + win32::GetLastErrorMessageA() );
#else
                                  + strerror( errno ) );
#endif
}

void write_rls_rle_level_set_file( const frantic::tstring& file, const levelset::rle_level_set& levelSet ) {
    if( levelSet.has_channel( _T("SignedDistance") ) )
        throw runtime_error( "write_rls_rle_level_set_file() - The level set being written to file \"" +
                             frantic::strings::to_string( file ) +
                             "\" has a \"SignedDistance\" channel, which is a reserved name." );

    frantic::tstring tempFilename = file + _T(".partial");

    FF_LOG( debug ) << "Writing parital file to " << tempFilename << endl;

    // This object will delete the temp file if an exception gets thrown.
    detail::rls_writer_tempfile_deleter tempFilenameDeleter( tempFilename );

    // Put the ofstream in a nested block
    {
        ofstream fout( tempFilename.c_str(), ios::out | ios::binary );

        if( !fout.is_open() )
            throw std::runtime_error( "write_rls_rle_level_set_file() - Failed to open the temporary file \"" +
                                      frantic::strings::to_string( tempFilename ) + "\" for output." );

        ////////
        // Write the main file header
        ////////

        // pass in an empty, fully inclusive channel policy
        detail::write_rle_header( fout, file, levelSet, channels::channel_propagation_policy(), true,
                                  levelSet.get_interface_voxel_width_inside(),
                                  levelSet.get_interface_voxel_width_outside() );

        ////////
        // Write the rle_index_spec
        ////////

        detail::write_rle_index_spec( fout, file, levelSet.get_rle_index_spec() );

        ////////
        // Write all the channel data
        ////////
        detail::write_rle_named_channels( fout, file, levelSet );

        ////////
        // Close the file, and make sure that worked
        ////////

        fout.close();
        if( fout.bad() )
            throw runtime_error( "write_rls_rle_level_set_file() - Failed to close the output stream for file \"" +
                                 frantic::strings::to_string( file ) + "\"." );
    }

    // Finally, Move the file we've been writing to its proper destination path
    // DeleteFile(file.c_str());
    // if( !MoveFile(tempFilename.c_str(), file.c_str()) )
    //	throw std::runtime_error( "write_rls_rle_level_set_file() - Failed to move the temporary file \"" + tempFilename
    //+
    //"\" to destination \"" + file + "\": " + win32::GetLastErrorMessage() );
    frantic::files::delete_file( file.c_str() );
    if( 0 != frantic::files::rename_file( tempFilename.c_str(), file.c_str() ) )
        throw std::runtime_error( "write_rle_voxel_field_file() - Failed to move the temporary file \"" +
                                  frantic::strings::to_string( tempFilename ) + "\" to destination \"" +
                                  frantic::strings::to_string( file ) + "\": "
#ifdef _WIN32
                                  + win32::GetLastErrorMessageA() );
#else
                                  + strerror( errno ) );
#endif
}

void interpolate_rls_rle_level_set_files( const frantic::tstring& firstFile, const frantic::tstring& secondFile,
                                          const float alpha, levelset::rle_level_set& outLevelSet ) {
    rle_level_set firstLevelSet, secondLevelSet;
    if( firstFile == secondFile )
        read_rle_level_set_file( firstFile, outLevelSet );
    else if( !( alpha > 0.f && alpha < 1.f ) )
        throw std::runtime_error( "interpolate_rls_rle_level_set_files() - Alpha value not between 0.0 and 1.0. (" +
                                  boost::lexical_cast<std::string>( alpha ) + ")" );
    else {
        read_rle_level_set_file( firstFile, firstLevelSet );
        read_rle_level_set_file( secondFile, secondLevelSet );
        outLevelSet.linear_interpolate( firstLevelSet, secondLevelSet, alpha );
    }
}

//
//
// rls_network_cache: DEPRECATED
//
//

// DEPRECATED
rls_network_cache::rls_network_cache( const frantic::tstring& tempDir, const frantic::tstring& cacheDir ) {
    initialize( tempDir, cacheDir );
}

// DEPRECATED
void rls_network_cache::initialize( const frantic::tstring& tempDir, const frantic::tstring& cacheDir ) {
    // init member variables
    m_tempDir = tempDir;
    m_fsq = frantic::files::filename_sequence( frantic::strings::to_tstring( cacheDir ) );
    m_fsq.sync_frame_set();
    m_initialized = true;

    // create the cache directory if it doesn't exist
    // if( !CreateDirectory( m_tempDir.c_str(), NULL ) ) { //null security attributes means it inherits the parent's
    //	if( GetLastError() != ERROR_ALREADY_EXISTS )
    //		throw runtime_error("rls_network_cache::initialize - Could not create temp directory for local caching ("
    //+ m_tempDir + ")");
    // }

    boost::filesystem::create_directory( boost::filesystem::path( m_tempDir ) );
}

// DEPRECATED
rls_network_cache::~rls_network_cache() {
    // remove all files copied by this object
    for( size_t i = 0; i < filesCopiedLocally.size(); ++i ) {
        const frantic::tstring& delFile = filesCopiedLocally[i];
        if( 0 != frantic::files::remove_file( delFile.c_str() ) ) {
            FF_LOG( warning ) << "rls_network_cache::~rls_network_cache() - Could not delete file: "
                              << frantic::strings::to_tstring( delFile ) << std::endl;
        }
    }
}

// DEPRECATED
bool rls_network_cache::get_level_set( double frame, float interfaceVoxelWidthInside, float interfaceVoxelWidthOutside,
                                       frantic::volumetrics::levelset::rle_level_set& outLevelSet,
                                       bool useCacheSettings ) {
    if( m_initialized ) {
        std::pair<double, double> bracket;
        float alpha;

        // grab the bracketing subframes and load the appropriate .rls files
        if( !get_nearest_subframe_interval( m_fsq, frame, bracket, alpha ) )
            throw runtime_error(
                "rls_network_cache::get_level_set - Could not locate an appropriate frame or frame bracket "
                "to load from cache " +
                frantic::strings::to_string( m_fsq.get_filename_pattern().get_pattern() ) + " for frame " +
                boost::lexical_cast<std::string, double>( frame ) );

        frantic::volumetrics::levelset::rle_level_set previousLevelSet, nextLevelSet;

        // if either of these files doesn't currently exist in our temp directory, copy them over
        // and check to make sure that they match the current scene level set (and each other) in terms
        // of coord system and interface widths
        frantic::tstring filename1, filename2, origFilename1, origFilename2;
        filename1 = m_tempDir + frantic::files::filename_from_path( ( m_fsq )[bracket.first] );

        origFilename1 = ( m_fsq )[bracket.first];

        if( !frantic::files::file_exists( filename1 ) ) {

            // psCopyCacheLevelSet.enter();
            copy_level_set_file( origFilename1, filename1, false );
            // psCopyCacheLevelSet.exit();

            read_rle_level_set_file_header( filename1, previousLevelSet );

            if( useCacheSettings ) {
                outLevelSet = volumetrics::levelset::rle_level_set( previousLevelSet.get_voxel_coord_system() );
            }

            if( !previousLevelSet.get_voxel_coord_system().equals( outLevelSet.get_voxel_coord_system(), 1e-5f ) ) {
                throw std::runtime_error(
                    "rls_network_cache::get_level_set - The voxel coordinate systems of the cached level "
                    "set and the scene are different(file: " +
                    frantic::strings::to_string( ( m_fsq )[bracket.first] ) +
                    " vcs: " + previousLevelSet.get_voxel_coord_system().str() +
                    ") versus (scene vcs: " + outLevelSet.get_voxel_coord_system().str() + ")." );
            }

            if( !useCacheSettings &&
                ( interfaceVoxelWidthInside > previousLevelSet.get_interface_voxel_width_inside() ||
                  interfaceVoxelWidthOutside > previousLevelSet.get_interface_voxel_width_outside() ) ) {
                std::string msg(
                    "rls_network_cache::get_level_set - The provided interface widths are greater than those "
                    "available in the cached level set (file " +
                    frantic::strings::to_string( ( m_fsq )[bracket.first] ) + " inside/outside interface widths: " +
                    boost::lexical_cast<std::string>( previousLevelSet.get_interface_voxel_width_inside() ) + "," +
                    boost::lexical_cast<std::string>( previousLevelSet.get_interface_voxel_width_outside() ) +
                    ") versus (scene: " + boost::lexical_cast<std::string>( interfaceVoxelWidthInside ) + "," +
                    boost::lexical_cast<std::string>( interfaceVoxelWidthOutside ) + ")." );
                cerr << "WARNING: " << msg << endl;
                cerr << "Not loading from cache" << endl;
                return false;
            }
        }

        if( !( alpha < 0 ) ) {
            filename2 = m_tempDir + frantic::files::filename_from_path( ( m_fsq )[bracket.second] );
            origFilename2 = ( m_fsq )[bracket.second];

            if( !frantic::files::file_exists( filename2 ) ) {
                // psCopyCacheLevelSet.enter();
                copy_level_set_file( origFilename2, filename2, false );
                // psCopyCacheLevelSet.exit();
                read_rle_level_set_file_header( filename2, nextLevelSet );

                if( useCacheSettings ) {
                    outLevelSet = volumetrics::levelset::rle_level_set( nextLevelSet.get_voxel_coord_system() );
                }

                if( !nextLevelSet.get_voxel_coord_system().equals( outLevelSet.get_voxel_coord_system(), 1e-5f ) )
                    throw std::runtime_error(
                        "rls_network_cache::get_level_set - The voxel coordinate systems of the cached "
                        "level set and the scene are different(file: " +
                        frantic::strings::to_string( ( m_fsq )[bracket.second] ) +
                        " vcs: " + nextLevelSet.get_voxel_coord_system().str() +
                        ") versus (scene vcs: " + outLevelSet.get_voxel_coord_system().str() + ")." );

                if( !useCacheSettings &&
                    ( interfaceVoxelWidthInside > nextLevelSet.get_interface_voxel_width_inside() ||
                      interfaceVoxelWidthOutside > nextLevelSet.get_interface_voxel_width_outside() ) ) {
                    std::string msg(
                        "rls_network_cache::get_level_set - The provided interface widths are greater than those "
                        "available in the cached level set (file " +
                        frantic::strings::to_string( ( m_fsq )[bracket.first] ) + " inside/outside interface widths: " +
                        boost::lexical_cast<std::string>( nextLevelSet.get_interface_voxel_width_inside() ) + "," +
                        boost::lexical_cast<std::string>( nextLevelSet.get_interface_voxel_width_outside() ) +
                        ") versus (scene: " + boost::lexical_cast<std::string>( interfaceVoxelWidthInside ) + "," +
                        boost::lexical_cast<std::string>( interfaceVoxelWidthOutside ) + ")." );
                    cerr << "WARNING: " << msg << endl;
                    cerr << "Not loading from cache" << endl;
                    return false;
                }
            }
        }

        if( alpha < 0 ) {
            read_rle_level_set_file( filename1, outLevelSet );
        } else {
            // psLoadingCacheLevelSet.enter();
            levelset::read_rle_level_set_file( filename1, previousLevelSet );
            // psLoadingCacheLevelSet.exit();
            // psLoadingCacheLevelSet.enter();
            levelset::read_rle_level_set_file( filename2, nextLevelSet );
            // psLoadingCacheLevelSet.exit();

            // interpolate the neighbouring level sets to create one at the current time
            // psLerpCacheLevelSet.enter();
            outLevelSet.linear_interpolate( previousLevelSet, nextLevelSet, alpha );
            // psLerpCacheLevelSet.exit();
        }

        // dump_profiling( cerr );
        return true;
    }

    return false; // not initialized
}

// DEPRECATED
bool rls_network_cache::get_nearest_subframe_interval( frantic::files::filename_sequence fsq, double frame,
                                                       std::pair<double, double>& interval, float& alpha ) {
    // Check if the file exists first
    if( frantic::files::file_exists( fsq[frame] ) ) {
        interval.first = frame;
        interval.second = frame;
        alpha = -1.0f; // this will serve as a flag stating no interpolation necessary
        return true;
    }

    // Grab the framenumbers from the cache
    vector<double> frameNumbers;
    fsq.allframe_numbers( frameNumbers );
    bool foundBefore = false, foundAfter = false;

    // Find the closest number before and after the given number in the filename sequence.
    double beforeNumber = boost::numeric::bounds<double>::lowest();
    double afterNumber = boost::numeric::bounds<double>::highest();

    for( vector<double>::iterator it = frameNumbers.begin(); it < frameNumbers.end(); ++it ) {

        if( frame - *it > 0 ) {
            if( frame - *it < frame - beforeNumber ) {
                beforeNumber = *it;
                foundBefore = true;
            }
        } else {
            if( *it - frame < afterNumber - frame ) {
                foundAfter = true;
                afterNumber = *it;
            }
        }
    }

    alpha = float( ( frame - beforeNumber ) / ( afterNumber - beforeNumber ) );

    interval.first = beforeNumber;
    interval.second = afterNumber;

    return ( foundBefore && foundAfter );
}

// DEPRECATED
void rls_network_cache::copy_level_set_file( const frantic::tstring& sourceFile, const frantic::tstring& destFile,
                                             const bool failFlag ) {
#ifdef _WIN32
    // try to copy the requested file
    if( CopyFile( sourceFile.c_str(), destFile.c_str(), failFlag ) )
        filesCopiedLocally.push_back( destFile );
    else
        throw runtime_error( "rls_network_cache::copy_level_set_file - Could not copy file to local temporary cache (" +
                             frantic::strings::to_string( sourceFile ) + " to " +
                             frantic::strings::to_string( destFile ) + ")" );
#endif
}

} // namespace levelset
} // namespace volumetrics
} // namespace frantic
