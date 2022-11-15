// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/files/aln.hpp>

#include <frantic/files/paths.hpp>
#include <frantic/graphics/transform4f.hpp>
#include <frantic/logging/logging_level.hpp>
#include <frantic/particles/particle_file_stream_factory.hpp>
#include <frantic/particles/streams/concatenated_particle_istream.hpp>
#include <frantic/particles/streams/set_channel_particle_istream.hpp>
#include <frantic/particles/streams/transformed_particle_istream.hpp>
#include <frantic/strings/tstring.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <fstream>
#include <iostream>
#include <utility>

using namespace frantic;
using namespace frantic::particles::streams;
using boost::format;
using boost::lexical_cast;
using boost::str;
using boost::filesystem::path;

namespace frantic {
namespace files {
namespace aln {

template <typename MatrixT>
void get_transforms( const tstring& file, std::vector<std::pair<tstring, MatrixT>>& outTransforms ) {
    std::ifstream fin( file.c_str() );

    if( fin.fail() )
        throw std::runtime_error( "aln.get_transforms: The file " + strings::to_string( file ) +
                                  " failed to open for reading." );

    outTransforms.clear();
    std::string scanName;
    MatrixT transform;
    int rowsValid = 0;  // How many rows in `transform` have been written?
    int lineNumber = 1; // Used for error messages
    int numScans = -1;
    // Define an "empty line" as a line with zero non-whitespace characters.
    // The outer loop upholds the following invariants:
    // 1. If `numScans == -1`, we are expecting to get the number of scans next or an empty line.
    // 2. If `scanName == ""`, we are expecting to get a file name, the EOF symbol "0", or an empty line next.
    // 3. If `scanName != "" && rowsValid == 0`, we are expecting either the "#" separator, or a row in a matrix next.
    // 4. If `scanName != "" && rowsValid != 0`, we are expecting a row in a matrix next.
    // If any of the above expectations are not met, you can expect an std::runtime_error to be thrown.
    for( std::string line; std::getline( fin, line ); ++lineNumber ) {
        boost::trim( line );
        if( line != "" ) {
            if( numScans == -1 ) {
                try {
                    numScans =
                        (int)lexical_cast<std::size_t>( line ); // convert to std::size_t first to catch negatives
                } catch( const boost::bad_lexical_cast& ) {
                    throw std::runtime_error( str(
                        format(
                            "aln.get_transforms: Found invalid number of scans \"%1%\" on line %2% of file \"%3%\"." ) %
                        line % lineNumber % strings::to_string( file ) ) );
                }
            } else if( scanName == "" ) {
                if( line != "0" ) {
                    // A file name
                    scanName = line;
                } else {
                    // End of file
                    break;
                }
            } else if( scanName != "" && rowsValid == 0 && line == "#" ) {
                // A separator character
                continue;
            } else {
                // Part of a matrix
                std::vector<std::string> row;
                boost::split( row, line, boost::is_any_of( " \t\v\f\r" ), boost::token_compress_on );

                if( row.size() != 4 ) {
                    throw std::runtime_error( str( format( "aln.get_transforms: Found an invalid number of transform "
                                                           "components on line %1% of file \"%2%\"." ) %
                                                   lineNumber % strings::to_string( file ) ) );
                }

                for( int i = 0; i < 4; ++i ) {
                    try {
                        // transform4t is column-major
                        transform[i * 4 + rowsValid] = lexical_cast<typename MatrixT::float_type>( row[i] );
                    } catch( const boost::bad_lexical_cast& ) {
                        throw std::runtime_error( str( format( "aln.get_transforms: Found invalid transform component "
                                                               "\"%1%\" on line %2% of file \"%3%\"." ) %
                                                       row[i] % lineNumber % strings::to_string( file ) ) );
                    }
                }
                rowsValid++;
                if( rowsValid == 4 ) {
                    outTransforms.push_back( std::make_pair( strings::to_tstring( scanName ), transform ) );
                    transform.set_to_zero();
                    scanName = "";
                    rowsValid = 0;
                }
            }
        }
    }

    if( numScans != outTransforms.size() ) {
        throw std::runtime_error( "aln.get_transforms: There number of scans specified in the header did not match the "
                                  "number of scans in the file body." );
    }
}

// Explicit template function instantiations
template void
get_transforms<graphics::transform4fd>( const tstring&,
                                        std::vector<std::pair<tstring, graphics::transform4fd>>& transforms );
template void
get_transforms<graphics::transform4f>( const tstring&,
                                       std::vector<std::pair<tstring, graphics::transform4f>>& transforms );

boost::shared_ptr<particle_istream> create_aln_particle_istream( const boost::filesystem::path& file,
                                                                 particles::particle_file_metadata& outMetadata,
                                                                 const channels::data_type_t positionTypeHint ) {
    std::vector<std::pair<tstring, graphics::transform4fd>> transforms;
    get_transforms<graphics::transform4fd>( strings::to_tstring( file.native() ), transforms );
    std::vector<tstring> scanNames;
    std::vector<graphics::transform4fd> matrices;
    for( size_t i = 0; i < transforms.size(); ++i ) {
        scanNames.push_back( transforms[i].first );
        matrices.push_back( transforms[i].second );
    }

    // Mostly copied from zfs_zfprj_utils.cpp
    std::vector<boost::shared_ptr<particle_istream>> scanStreams;
    frantic::particles::particle_file_stream_factory_object factory;
    factory.set_position_type_hint( positionTypeHint );
    for( std::size_t i = 0; i < scanNames.size(); ++i ) {
        particle_istream_ptr rawStream =
            ( path( scanNames[i] ).is_absolute()
                  ? factory.create_istream( scanNames[i] )
                  : factory.create_istream( files::to_tstring( file.parent_path() / scanNames[i] ) ) );
        particle_istream_ptr transformedStream( new transformed_particle_istream<double>( rawStream, matrices[i] ) );
        scanStreams.push_back( boost::make_shared<set_channel_particle_istream<boost::uint32_t>>(
            transformedStream, _T("ScannerIndex"), static_cast<boost::uint32_t>( i + 1 ) ) );
    }

    channels::property_map generalMetadata;
    channels::channel_map generalMap;
    particles::prt::add_coordinate_system( generalMap );
    particles::prt::length_unit_in_micrometers::add_channel( generalMap );
    generalMap.end_channel_definition();
    generalMetadata.set_channel_map( generalMap );
    particles::prt::set_coordinate_system( generalMetadata, frantic::graphics::coordinate_system::right_handed_zup );
    particles::prt::length_unit_in_micrometers::set_value( generalMetadata, 1e6 );
    particles::prt::set_scanner_transforms( generalMetadata, matrices );
    outMetadata.set_general_metadata( generalMetadata );

    boost::shared_ptr<particle_istream> concatenatedStream( new concatenated_particle_istream( scanStreams ) );
    // Set to native channel map so we get the "ScannerIndex" channel.
    concatenatedStream->set_channel_map( concatenatedStream->get_native_channel_map() );
    return concatenatedStream;
}
} // namespace aln
} // namespace files
} // namespace frantic
