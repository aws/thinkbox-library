// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
/**
 * @file
 * @author Stephen Kiazyk
 *
 * Potentially useful extensions to the <utility> header
 */

#pragma once

#include <algorithm>
#include <utility>

namespace frantic {

/**
 * Returns a homogeneous pair that is in 'sorted' order.
 * This is useful if you need a pair type as elements of a set.
 */
template <class T>
std::pair<T, T> make_sorted_pair( const T& u, const T& v ) {
    return std::make_pair( std::min( u, v ), std::max( u, v ) );
}

template <class T>
inline void clear_with_swap( std::vector<T>& v ) {
    std::vector<T> temp;
    v.swap( temp );
}

template <class T>
inline void shrink_to_fit( std::vector<T>& v ) {
    std::vector<T> temp( v.begin(), v.end() );
    v.swap( temp );
}

} // namespace frantic
