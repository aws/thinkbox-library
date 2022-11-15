// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace frantic {
namespace strings {

bool is_valid_utf8( const char* sbegin, const char* send );
bool is_valid_utf8( const char* s );
bool is_valid_utf8( const std::string& s );

std::wstring wstring_from_utf8( const char* s );
std::wstring wstring_from_utf8( const std::string& s );

std::string to_utf8( const char* s );
std::string to_utf8( const std::string& s );

std::string to_utf8( const wchar_t* s );
std::string to_utf8( const std::wstring& s );

} // namespace strings
} // namespace frantic
