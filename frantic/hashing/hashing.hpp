// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace frantic {
namespace hashing {

/**
 * \param out a pointer to hash_size bytes to be uses as output
 * \param in a pointer to the bytes you want to hash
 * \param inlen the number of bytes to hash
 */
template <class HashAlgorithm>
inline void get_hash( unsigned char* out, const void* in, std::size_t inlen ) {
    HashAlgorithm h;
    h.update( in, inlen );
    h.final_value( out );
}

template <class HashAlgorithm>
inline std::string get_hash( const void* in, std::size_t inlen ) {
    unsigned char* out = (unsigned char*)alloca( HashAlgorithm::hash_size );
    get_hash<HashAlgorithm>( out, in, inlen );
    return hash_to_string( out, HashAlgorithm::hash_size );
}

/**
 * Convert a hexadecimal hash to a string.
 *
 * @param in a pointer to the bytes to hash.
 * @param inSize the number of bytes in the hash.
 */
std::string hash_to_string( unsigned char* in, std::size_t inSize );

} // namespace hashing
} // namespace frantic
