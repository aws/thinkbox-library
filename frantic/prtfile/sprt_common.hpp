// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

// NOTE: This here is assuming a little-endian machine, as do the PRT and PRT2 I/O routines.

// {0xC0, 't', 's', 't', '\r', '\n', 0x1A, '\n'};
#define SPRT1_MAGIC_NUMBER ( 0x0A1A0A0D747374C0ULL )

namespace frantic {
namespace prtfile {

inline bool has_sprt1_magic_number( const frantic::tstring& filename ) {
    std::ifstream fin( filename.c_str(), std::ios::in | std::ios::binary );
    boost::uint64_t magic;
    fin.read( reinterpret_cast<char*>( &magic ), 8u );
    return ( !fin.fail() && magic == SPRT1_MAGIC_NUMBER );
}

inline std::size_t get_expected_number_of_octree_levels( boost::uint64_t totalCount, boost::uint64_t lowestDetailTotal,
                                                         std::size_t levelCountMultiplier ) {
    std::size_t numLevels;
    if( totalCount <= lowestDetailTotal ) {
        numLevels = 1;
    } else {
        numLevels = static_cast<std::size_t>( std::ceil(
            ( std::log( static_cast<double>( totalCount ) ) - std::log( static_cast<double>( lowestDetailTotal ) ) ) /
            std::log( static_cast<double>( levelCountMultiplier ) ) ) );
    }
    return numLevels;
}

} // namespace prtfile
} // namespace frantic
