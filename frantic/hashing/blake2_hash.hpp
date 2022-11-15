// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#pragma warning( push )
#pragma warning( disable : 4804 )
#include <libb2/blake2.h>
#pragma warning( pop )

namespace frantic {
namespace hashing {

/*
 * BLAKE2s is optimized for 8- to 32-bit platforms and produces digests
 * of any size between 1 and 32 bytes
 */
class blake2s_hash {
  public:
    enum { hash_size = BLAKE2S_OUTBYTES };

    blake2s_hash();
    void update( const void* in, boost::uint64_t inlen );
    void final_value( unsigned char* out );

  private:
    blake2s_state m_state;
};

/*
 * BLAKE2b (or just BLAKE2) is optimized for 64-bit platforms - including
 * NEON-enabled ARMs - and produces digests of any size between 1 and 64 bytes
 */
class blake2b_hash {
  public:
    enum { hash_size = BLAKE2B_OUTBYTES };

    blake2b_hash();
    void update( const void* in, boost::uint64_t inlen );
    void final_value( unsigned char* out );

  private:
    blake2b_state m_state;
};

/*
 * 8-way parallel BLAKE2sp
 */
class blake2sp_hash {
  public:
    enum { hash_size = BLAKE2S_OUTBYTES };

    blake2sp_hash();
    void update( const void* in, boost::uint64_t inlen );
    void final_value( unsigned char* out );

  private:
    blake2sp_state m_state;
};

/*
 * 4-way parallel BLAKE2bp
 */
class blake2bp_hash {
  public:
    enum { hash_size = BLAKE2B_OUTBYTES };

    blake2bp_hash();
    void update( const void* in, boost::uint64_t inlen );
    void final_value( unsigned char* out );

  private:
    blake2bp_state m_state;
};

} // namespace hashing
} // namespace frantic
