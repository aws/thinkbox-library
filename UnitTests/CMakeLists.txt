# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
cmake_minimum_required( VERSION 3.15 FATAL_ERROR )

project( ThinkboxLibraryTest )

find_package( thinkboxcmlibrary REQUIRED )
include( PrecompiledHeader)
include( ThinkboxCMLibrary)

add_executable( test_thinkboxlibrary "" )

set_property( TARGET test_thinkboxlibrary PROPERTY CXX_STANDARD 17 )

target_include_directories( test_thinkboxlibrary
	PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} )

# Treat warnings as errors
if( WIN32 )
	target_compile_options(test_thinkboxlibrary PRIVATE "/W3" "/WX")
endif()

target_include_directories( test_thinkboxlibrary PRIVATE "../include" )
target_link_libraries( test_thinkboxlibrary PRIVATE thinkboxlibrary )

find_package( GTest REQUIRED )
target_include_directories( test_thinkboxlibrary PRIVATE ${GTest_INCLUDE_DIRS} )
target_link_libraries( test_thinkboxlibrary PRIVATE GTest::GTest )

# The Conan version of Boost was built with this, and it changes the library names.
# As a result, we need to set this to tell Boost to look for the right libraries to
# link against.
# target_compile_definitions( test_thinkboxlibrary PUBLIC BOOST_AUTO_LINK_SYSTEM )

file( GLOB_RECURSE H_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "*.h"
      "*.hpp"
)

file( GLOB_RECURSE CXX_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "*.cpp"
      "*.cc"
)

target_sources( test_thinkboxlibrary PRIVATE
	stdafx.cpp
	stdafx.h
	${H_FILES}
	${CXX_FILES}
)

frantic_common_platform_setup( test_thinkboxlibrary )
frantic_default_source_groups( test_thinkboxlibrary )
frantic_link_apple_core_libraries( test_thinkboxlibrary )

add_precompiled_header( test_thinkboxlibrary stdafx.h SOURCE_CXX stdafx.cpp )

# Disable optimization for the RelWithDebInfo configuration on Windows.
# This allows breakpoints to be hit reliably when debugging in Visual Studio.
if( WIN32 )
	target_compile_options( test_thinkboxlibrary PRIVATE "$<$<CONFIG:Release>:/O2>$<$<CONFIG:RelWithDebInfo>:/Od>" )
endif()
