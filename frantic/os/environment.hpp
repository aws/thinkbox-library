// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/strings/tstring.hpp>

namespace frantic {
namespace os {

/**
 * Set the environment variable with name key to value for the duration of current process.
 *
 * @param key the name of the environment variable
 * @param value the value to store
 *
 * Will throw a runtime error if it failed to set the environment variable.
 */
void set_environment_variable( const frantic::tstring& key, const frantic::tstring& value );

/**
 * Return the environment variable with name key.
 *
 * @param key the name of the environment variable
 *
 * @return the value stored in the environment variable if it's set, an empty string otherwise
 */
frantic::tstring get_environment_variable( const frantic::tstring& key );

/**
 * Return the environment variable with name key.
 *
 * Unlike get_environment_variable, this will attempt to read the most up to date environment variable to account
 * for them changing while the program is running.  It is less efficient than get_environment_variable as it needs
 * to construct the updated list of variables and iterate through it to find the entry.
 *
 * This is only implemented for Windows.  On other OSs, it just calls get_environment_variable.
 *
 * If the environment block could not be initialized, it will call get_environment_variable.
 *
 * @param key the name of the environment variable
 * @param outIsFromOS if not NULL, this is set to true if it tried to read the value from the OS environment
 *                    variable block.  Otherwise it is set to false.
 *
 * @return the value stored in the environment variable if it's set, an empty string otherwise
 */
frantic::tstring get_environment_variable_from_os_scope( const frantic::tstring& key, bool* outIsFromOS = NULL );

} // namespace os
} // namespace frantic
