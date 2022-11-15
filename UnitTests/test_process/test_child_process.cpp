// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#include <frantic/os/child_process.hpp>

using std::string;

TEST( ChildProcess, EchoHelloWorld ) {
    frantic::process::child_process process;
    process.set_controlStdOut( true );
#ifdef _WIN32
    process.launch( _T( "cmd" ), _T( "/C echo Hello World" ) );
#else
    process.launch( "echo", "Hello World" );
#endif

    string output = "";
    EXPECT_TRUE( process.getstdoutline( -1, output ) );
    EXPECT_EQ( "Hello World", output );
    bool complete = process.waitforexit();
    ASSERT_TRUE( complete );
    int exitcode = process.getexitcode();
    EXPECT_EQ( 0, exitcode );
}
