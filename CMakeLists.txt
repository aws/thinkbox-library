# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

project( ThinkboxLibrary )

find_package( thinkboxcmlibrary REQUIRED )
include( PrecompiledHeader )
include( ThinkboxCMLibrary )

option( BUILD_UNIT_TESTS "Build unit tests" ON )

option( BUILD_WITH_TBB "Build with Intel Thread Building Blocks" ON )
option( BUILD_WITH_E57 "Build with xerces and libe57 to enable e57 file IO" ON )
option( BUILD_WITH_LZ4 "Build with LZ4 to enable PRT2 file format" ON )

add_library( thinkboxlibrary STATIC )

set_property( TARGET thinkboxlibrary PROPERTY CXX_STANDARD 17 )

target_include_directories(thinkboxlibrary PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>
)

file( GLOB_RECURSE H_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "frantic/*.h"
      "frantic/*.hpp"
)

file( GLOB_RECURSE CXX_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
      "src/*.cpp"
)

target_sources( thinkboxlibrary PRIVATE
	stdafx.cpp
	stdafx.h
    ${H_FILES}
    ${CXX_FILES}
)

# Treat warnings as errors
if(WIN32)
	target_compile_options(thinkboxlibrary PRIVATE "/W3" "/WX")
	# Disable warning C4503: '<symbol>' : decorated name length exceeded, name was truncated
	target_compile_options(thinkboxlibrary PUBLIC "/wd4503")
else()
	if(CMAKE_COMPILER_IS_GNUCXX)
		# Disable warning: unique_ptr is deprecated
		target_compile_options(thinkboxlibrary PUBLIC "-Wno-deprecated-declarations")
		# Disable warning: assuming signed overflow does not occur when assuming that (X + c) < X is always false
		target_compile_options(thinkboxlibrary PUBLIC "-Wno-strict-overflow")
		# Because of https://gcc.gnu.org/bugzilla/show_bug.cgi?id=53431, can't selectively remove this from boost
		target_compile_options(thinkboxlibrary PUBLIC "-Wno-unused-variable")
	endif()
endif()

find_package( libb2 REQUIRED )
find_package( Boost REQUIRED )
find_package( BZip2 REQUIRED )
find_package( Eigen3 REQUIRED )
find_package( glog REQUIRED )
find_package( OpenEXR REQUIRED )
find_package( ZLIB REQUIRED )
find_package( TBB REQUIRED )
find_package( tinyxml2 REQUIRED )
find_package( utf8cpp REQUIRED )
find_package( xxHash REQUIRED )

target_include_directories( thinkboxlibrary PUBLIC ${libb2_INCLUDE_DIRS} )
target_include_directories( thinkboxlibrary PUBLIC ${Boost_INCLUDE_DIRS} )
target_include_directories( thinkboxlibrary PUBLIC ${BZZip2_INCLUDE_DIRS} )
target_include_directories( thinkboxlibrary PUBLIC ${Eigen3_INCLUDE_DIRS} )
target_include_directories( thinkboxlibrary PUBLIC ${glog_INCLUDE_DIRS} )
target_include_directories( thinkboxlibrary PUBLIC ${OpenEXR_INCLUDE_DIRS} )
target_include_directories( thinkboxlibrary PUBLIC ${ZLIB_INCLUDE_DIRS} )
target_include_directories( thinkboxlibrary PUBLIC ${TBB_INCLUDE_DIRS} )
target_include_directories( thinkboxlibrary PUBLIC ${tinyxml2_INCLUDE_DIRS} )
target_include_directories( thinkboxlibrary PUBLIC ${utf8cpp_INCLUDE_DIRS} )
target_include_directories( thinkboxlibrary PUBLIC ${xxHash_INCLUDE_DIRS} )

target_link_libraries( thinkboxlibrary INTERFACE libb2::libb2 )
target_link_libraries( thinkboxlibrary INTERFACE Boost::Boost )
target_link_libraries( thinkboxlibrary INTERFACE BZip2::BZip2 )
target_link_libraries( thinkboxlibrary INTERFACE Eigen3::Eigen )
target_link_libraries( thinkboxlibrary INTERFACE glog::glog )
target_link_libraries( thinkboxlibrary INTERFACE OpenEXR::OpenEXR )
target_link_libraries( thinkboxlibrary INTERFACE ZLIB::ZLIB )
target_link_libraries( thinkboxlibrary INTERFACE TBB::tbb )
target_link_libraries( thinkboxlibrary INTERFACE tinyxml2::tinyxml2 )
target_link_libraries( thinkboxlibrary INTERFACE utf8cpp::utf8cpp )
target_link_libraries( thinkboxlibrary INTERFACE xxHash::xxHash )

if( UNIX AND NOT APPLE )
	find_package( ICU REQUIRED )
	target_include_directories( thinkboxlibrary PUBLIC ${ICU_INCLUDE_DIRS})
	target_link_libraries( thinkboxlibrary INTERFACE ICU::ICU )
endif()

if( BUILD_WITH_E57 )
	find_package( XercesC REQUIRED )
	find_package( e57format REQUIRED )
	target_include_directories( thinkboxlibrary PUBLIC ${XercesC_INCLUDE_DIRS} )
	target_include_directories( thinkboxlibrary PUBLIC ${E57Format_INCLUDE_DIRS} )
	target_link_libraries( thinkboxlibrary INTERFACE XercesC::XercesC )
	target_link_libraries( thinkboxlibrary INTERFACE E57Format::E57Format )
endif()
if( BUILD_WITH_LZ4 )
	find_package( lz4 REQUIRED )
	target_include_directories( thinkboxlibrary PUBLIC ${lz4_INCLUDE_DIRS} )
	target_link_libraries( thinkboxlibrary INTERFACE lz4::lz4 )
endif()


frantic_default_source_groups( thinkboxlibrary HEADERDIR thinkboxlibrary SOURCEDIR src )
frantic_common_platform_setup( thinkboxlibrary )

add_precompiled_header( thinkboxlibrary stdafx.h SOURCE_CXX stdafx.cpp )

# Disable optimization for the RelWithDebInfo configuration on Windows.
# This allows breakpoints to be hit reliably when debugging in Visual Studio.
if( WIN32 )
	target_compile_options( thinkboxlibrary PRIVATE "$<$<CONFIG:Release>:/O2>$<$<CONFIG:RelWithDebInfo>:/Od>" )
endif()

# The Conan version of Boost was built with this, and it changes the library names.
# As a result, we need to set this to tell Boost to look for the right libraries to
# link against.
target_compile_definitions( thinkboxlibrary PUBLIC BOOST_AUTO_LINK_SYSTEM )

# We have to use boost 1.78 for now because the version of openimageio in Conan isn't compatible with boost 1.79.
# In order to get it to compile we need to define the following so C++17 deprecated allocators don't cause a compiler
# warning which will be interpreted as an error.
target_compile_definitions( thinkboxlibrary PUBLIC _SILENCE_CXX17_OLD_ALLOCATOR_MEMBERS_DEPRECATION_WARNING )

# Cause a compiler error to be thrown if any of Eigen's LGPL licensed code is included to
# avoid accidentally using copy-left code.
target_compile_definitions( thinkboxlibrary PUBLIC EIGEN_MPL2_ONLY )

if( BUILD_WITH_E57 )
	target_compile_definitions( thinkboxlibrary PUBLIC E57_AVAILABLE )
endif()

if( BUILD_WITH_LZ4 )
	target_compile_definitions( thinkboxlibrary PUBLIC LZ4_AVAILABLE )
endif()

if( NOT BUILD_WITH_TBB )
	target_compile_definitions( thinkboxlibrary PUBLIC FRANTIC_DISABLE_THREADS )
else()
	target_compile_definitions( thinkboxlibrary PUBLIC FRANTIC_TBB_AVAILABLE )
endif()

if( BUILD_UNIT_TESTS )
	add_subdirectory( UnitTests )
endif()

install( DIRECTORY frantic
		 DESTINATION include
		 FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)
install( TARGETS thinkboxlibrary
		 RUNTIME DESTINATION bin
		 LIBRARY DESTINATION lib
		 ARCHIVE DESTINATION lib
)
