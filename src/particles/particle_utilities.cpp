// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/particles/particle_utilities.hpp>
#include <frantic/particles/streams/csv_particle_istream.hpp>

#include <boost/algorithm/string.hpp>

namespace frantic {
namespace particles {

bool is_csv_particle_file( const frantic::tstring& filename ) {

    particle_file_stream_factory_object factory;

    bool found = false;

    frantic::tstring ext = frantic::files::extension_from_path( filename );
    boost::to_lower( ext );

    // Method 1: check the extension (presumably this could be extended to other known csv format extensions)
    found = ( ext == _T(".csv") );

    if( !found ) {

        try {
            // Method 2: check if the stream returned is a csv_particle_istream
            particle_istream_ptr p = factory.create_istream( filename );
            found = dynamic_cast<streams::csv_particle_istream*>( p.get() ) != NULL;
        } catch( std::runtime_error& ) {
            found = false;
        }
    }

    return found;
}

bool is_binary_particle_file( const frantic::tstring& filename ) {
    frantic::tstring ext = frantic::files::extension_from_path( filename );
    boost::to_lower( ext );

    bool isBinary = ext == _T( ".prt" ) || ext == _T( ".sprt" ) ||
#if defined( E57_AVAILABLE )
                    ext == _T( ".e57" ) ||
#endif
                    ext == _T( ".las" ) || ext == _T( ".ptg" ) || ext == _T( ".bin" ) || ext == _T( ".rpc" );

    return isBinary;
}

} // namespace particles
} // namespace frantic
