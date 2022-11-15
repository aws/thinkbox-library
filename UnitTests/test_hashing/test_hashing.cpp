// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <boost/cstdint.hpp>

#include <tbb/task_scheduler_init.h>

#include <frantic/hashing/blake2_hash.hpp>
#include <frantic/hashing/hashing.hpp>

using namespace frantic::hashing;

TEST( Hashing, HashToString ) {
    boost::uint8_t zero[1] = { 0 };
    EXPECT_EQ( "00", hash_to_string( zero, sizeof( zero ) ) );

    boost::uint8_t abcd[2] = { 0xAB, 0xCD };
    EXPECT_EQ( "abcd", hash_to_string( abcd, sizeof( abcd ) ) );
}

TEST( Hashing, Blake2b ) {
    // example from RFC 7693, "BLAKE2 Crypto Hash and MAC"
    const unsigned char abc[] = { 'a', 'b', 'c' };
    const std::string expected = "ba80a53f981c4d0d6a2797b69f12f6e9"
                                 "4c212f14685ac4b74b12bb6fdbffa2d1"
                                 "7d87c5392aab792dc252d5de4533cc95"
                                 "18d38aa8dbf1925ab92386edd4009923";
    EXPECT_EQ( expected, get_hash<blake2b_hash>( abc, sizeof( abc ) ) );
}

TEST( Hashing, Blake2bp ) {
    // hash I generated myself, so it shouldn't necessarily be trusted
    const unsigned char abc[] = { 'a', 'b', 'c' };
    const std::string expected = "b91a6b66ae87526c400b0a8b53774dc6"
                                 "5284ad8f6575f8148ff93dff943a6ecd"
                                 "8362130f22d6dae633aa0f91df4ac89a"
                                 "aff31d0f1b923c898e82025dedbdad6e";
    EXPECT_EQ( expected, get_hash<blake2bp_hash>( abc, sizeof( abc ) ) );
}

TEST( Hashing, Blake2bpLarge ) {
    // Hash an array large enough to use the multithreaded code path
    tbb::task_scheduler_init taskScheduler;

    // hash I generated myself, so it shouldn't necessarily be trusted
    std::vector<unsigned char> zeros( 1000000 );
    const std::string expected = "b56224e79b8305fc7b2045ef9fd02f4d"
                                 "1ed97e8b170fb409d03e12d28691b23e"
                                 "08952e34539c3265c8f98251118bca91"
                                 "c664d12924610a77400958772f2ca579";
    EXPECT_EQ( expected, get_hash<blake2bp_hash>( &zeros[0], zeros.size() ) );
}

TEST( Hashing, Blake2s ) {
    // example from RFC 7693, "BLAKE2 Crypto Hash and MAC"
    const unsigned char abc[] = { 'a', 'b', 'c' };
    const std::string expected = "508c5e8c327c14e2e1a72ba34eeb452f"
                                 "37458b209ed63a294d999b4c86675982";
    EXPECT_EQ( expected, get_hash<blake2s_hash>( abc, sizeof( abc ) ) );
}

TEST( Hashing, Blake2sp ) {
    // hash I generated myself, so it shouldn't necessarily be trusted
    const unsigned char abc[] = { 'a', 'b', 'c' };
    const std::string expected = "70f75b58f1fecab821db43c88ad84edd"
                                 "e5a52600616cd22517b7bb14d440a7d5";
    EXPECT_EQ( expected, get_hash<blake2sp_hash>( abc, sizeof( abc ) ) );
}

TEST( Hashing, Blake2spLarge ) {
    // Hash an array large enough to use the multithreaded code path
    tbb::task_scheduler_init taskScheduler;

    // hash I generated myself, so it shouldn't necessarily be trusted
    std::vector<unsigned char> zeros( 1000000 );
    const std::string expected = "175ce84373591fdd19a9eeec7fd7e3ae"
                                 "a74eb3b1ee5d42d94a9ce6218c315f52";
    EXPECT_EQ( expected, get_hash<blake2sp_hash>( &zeros[0], zeros.size() ) );
}
