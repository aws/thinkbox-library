// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
// clang-format off
#include "stdafx.h"
// clang-format on

#if defined( _WIN32 )

#include <boost/algorithm/string.hpp>
#include <boost/uuid/string_generator.hpp>

#include "gtest/gtest.h"

#include <frantic/os/child_process.hpp>
#include <frantic/win32/utility.hpp>

TEST( WIN_32, GetUUID ) {
    /**
     * Tests that get_win32_computer_system_uuid returns the same UUID as 'wmic path win32_computersystemproduct get
     * uuid'
     */
    frantic::process::child_process cp;
    cp.set_terminateOnExit( true );
    cp.set_controlStdOut( true );
    cp.set_hideWindow( true );

    /**
     * wmic is installed on all windows operating systems
     * https://docs.microsoft.com/en-us/windows/win32/wmisdk/operating-system-availability-of-wmi-components
     * With win32_computersystemproduct being on operating systems Windows Vista and higher
     * https://docs.microsoft.com/en-us/windows/win32/cimwin32prov/win32-computersystemproduct
     */
    cp.launch( _T("cmd"), _T("/C wmic path win32_computersystemproduct get uuid"), _T("./") );

    std::string output;
    // first line of command output is "UUID"
    EXPECT_TRUE( cp.getstdoutline( 1000, output ) );
    // second line of command output is a blank line
    EXPECT_TRUE( cp.getstdoutline( 1000, output ) );
    // third line of command output is the actual UUID value such as "EC2AE145-D1DC-13B2-94ED-01234ABCDEF    "
    EXPECT_TRUE( cp.getstdoutline( 1000, output ) );

    // the UUID has spaces on the end, trim them off
    boost::trim_right( output );

    frantic::win32::get_win32_computer_system_uuid get_uuid;
    const std::wstring uuid = get_uuid.get();
    boost::uuids::string_generator gen;
    EXPECT_EQ( gen( uuid ), gen( output ) );
    EXPECT_EQ( gen( uuid ),
               gen( output ) ); // calling again to ensure CoInitializeEx and CoUninitialize were called correctly
}
#endif