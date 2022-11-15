// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/hashing/hashing.hpp>

namespace frantic {
namespace hashing {

namespace {

char get_hex_char_from_nybble( unsigned char nybble ) {
    assert( nybble < 16 );
    const char lut[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };
    return lut[nybble];
}

} // anonymous namespace

std::string hash_to_string( unsigned char* in, std::size_t inSize ) {
    std::string result;
    result.resize( 2 * inSize );

    for( std::size_t i = 0; i < inSize; i++ ) {
        const unsigned char c = in[i];
        result[2 * i] = get_hex_char_from_nybble( c >> 4 );
        result[2 * i + 1] = get_hex_char_from_nybble( c & 0x0F );
    }

    return result;
}

} // namespace hashing
} // namespace frantic
