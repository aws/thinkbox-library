// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <stdexcept>

// macros for nice error messages
// if the assert condition fails,
// FRANTIC_ASSERT_THROW() will throw a standard runtime error
// with the message "FileName:LineNumber FunctionName(): assert condition " +
//      a string representation of the assert condition + " failed: " + the message you supply.
namespace frantic {
namespace diagnostics {
namespace detail {

inline void assert_throw_runtime_error( bool assert_condition, std::string assert_condition_string,
                                        std::string failure_message, std::string file_name, std::string line_number,
                                        std::string function_name ) {
    if( !assert_condition ) {
        std::string err_message = file_name.substr( file_name.find_last_of( '\\' ) + 1 ) + ":" + line_number + " " +
                                  function_name + "() assert condition " + assert_condition_string +
                                  " failed: " + failure_message;
        throw std::runtime_error( err_message );
    }
}

} // namespace detail
} // namespace diagnostics
} // namespace frantic

#define FRANTIC_INTERNAL_STRINGIFY( x ) #x
#define FRANTIC_INTERNAL_TOSTRING( x ) FRANTIC_INTERNAL_STRINGIFY( x )

#define FRANTIC_ASSERT_THROW( assert_condition, failure_message )                                                      \
    frantic::diagnostics::detail::assert_throw_runtime_error( assert_condition, #assert_condition, failure_message,    \
                                                              __FILE__, FRANTIC_INTERNAL_TOSTRING( __LINE__ ),         \
                                                              __FUNCTION__ )

// a debug version for asserts not wanted in the release version
#ifdef NDEBUG
#define FRANTIC_DEBUG_ASSERT_THROW( assert_condition, failure_message )
#else
#define FRANTIC_DEBUG_ASSERT_THROW( assert_condition, failure_message )                                                \
    frantic::diagnostics::detail::assert_throw_runtime_error( assert_condition, #assert_condition, failure_message,    \
                                                              __FILE__, FRANTIC_INTERNAL_TOSTRING( __LINE__ ),         \
                                                              __FUNCTION__ )
#endif
