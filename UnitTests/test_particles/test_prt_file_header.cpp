// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include "gtest/gtest.h"

#include <fstream>

#include <frantic/misc/exception_stream.hpp>
#include <frantic/particles/prt_file_header.hpp>

// To dump the input files for these tests you can use a utility called od (included with git bash).
// `od -t x1 <input_filename>` can be used to dump the contents as hexadecimal.
// `od -t c <input_filename>` can be used to dump the contents as printable characters or backslash escape.
// See `od --help` for a full list of options.

TEST( PRTFileHeader, FakeMassiveChunkDivisibleInfiniteLoop ) {
    // If a PRT file specifies a chunk size that would result in the end of the chunk being beyond the end of the file
    // and then you use seekg, this can cause the boolean operator for the stream to return true even if
    // the stream is technically at EOF (this seems to be related to the implementation, as I was unable to reproduce
    // using vc14). This is a regression test to ensure that this infinite looping behavior does not get re-introduced
    // for a PRT file with specifications matching what we were given in a pentest finding.
    frantic::particles::prt_file_header header;
    frantic::tstring streamName( _T( "TestInputs/massive_chunk_divisible.prt" ) );
    std::ifstream inputStream( "TestInputs/massive_chunk_divisible.prt", std::ios::binary );

    // If the test fails, this call will cause an infinite loop.
    // It will timeout in the context of CI, and will otherwise need to
    // manually be stopped.
    EXPECT_THROW( header.read_header( inputStream, streamName ), std::runtime_error );
}

TEST( PRTFileHeader, FakeMassiveChunkIndivisibleInfiniteLoop ) {
    // Similar to the above test, however, in this case the behaviour is caught earlier due to a check to make sure the
    // length of a chunk data section is divisible by the chunk data type's size.
    frantic::particles::prt_file_header header;
    frantic::tstring streamName( _T( "TestInputs/massive_chunk_indivisible.prt" ) );
    std::ifstream inputStream( "TestInputs/massive_chunk_indivisible.prt", std::ios::binary );

    // If the test fails, this call will cause an infinite loop.
    // It will timeout in the context of CI, and will otherwise need to
    // manually be stopped.
    EXPECT_THROW( header.read_header( inputStream, streamName ), std::runtime_error );
}

TEST( PRTFileHeader, NegativeEightChunksizeInfiniteLoop ) {
    // If the reader doesn't recognize the type of a chunk, it will skip that chunk.
    // To skip the chunk it will move forward the size of the chunk.
    // If the chunk size is set to -8 this will cause an infinite loop by trying to read the same chunk
    // over and over again.
    frantic::particles::prt_file_header header;
    frantic::tstring streamName( _T( "TestInputs/negative_eight.prt" ) );
    std::ifstream inputStream( "TestInputs/negative_eight.prt", std::ios::binary );

    // If the test fails, this call will cause an infinite loop.
    // It will timeout in the context of CI, and will otherwise need to
    // manually be stopped.
    EXPECT_THROW( header.read_header( inputStream, streamName ), std::runtime_error );
}
